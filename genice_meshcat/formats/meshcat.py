# coding: utf-8

desc={"ref": {"Meshcat": "https://github.com/rdeits/meshcat-python"},
      "brief": "Meshcat + Google Colab.",
      "usage": ""}


from collections import defaultdict
from logging import getLogger
import numpy as np
from genice2.decorators import timeit, banner
import genice2.formats
import meshcat
from meshcat.jupyter import JupyterVisualizer
import meshcat.geometry as g
import meshcat.transformations as tf
import numpy as np
from math import acos, pi
import yaplotlib as yp # transient
from time import sleep


def draw_atom(v, label, atom, radius, color=0xffffff):
    v[label].set_object(g.Mesh(g.Sphere(radius), g.MeshLambertMaterial(color=color)))
    v[label].set_transform(tf.translation_matrix(atom))


def draw_bond(v, label, p1, d, radius, color=0xffffff):
    H = np.linalg.norm(d)
    R = np.linalg.norm(d[:2])
    e = d/H
    x = np.array([0,1,0])
    rot  = -acos(e@x)
    if -1 < e@x < 1:
        axis = np.cross(e,x)
        v[label].set_object(g.Mesh(g.Cylinder(H, radius), g.MeshLambertMaterial(color=color)))
        v[label].set_transform(tf.translation_matrix(p1).dot(tf.rotation_matrix(rot, axis).dot(tf.translation_matrix([0, H/2, 0]))))
    else:
        v[label].set_object(g.Mesh(g.Cylinder(H, radius), g.MeshLambertMaterial(color=color)))
        v[label].set_transform(tf.translation_matrix(p1+d/2))



class Format(genice2.formats.Format):
    """
Output the atomic positions on Colab using Meshcat.

Options:
    H=x   Set the radius of H to be x
    """


    def __init__(self, **kwargs):
        unknown = dict()
        self.size_H = 0.015
        jupyter = False
        for k, v in kwargs.items():
            if k == "H":
                self.size_H = float(v)
            elif k == "Jupyter":
                if v is None:
                    jupyter = True
                elif v is True:
                    jupyter = True
            else:
                unknown[k] = v
        super().__init__(**unknown)
        self.output = None
        if Jupyter:
            # Jupyter and colab inline
            self.vis = JupyterVisualizer()
        else:
            # On a Browser window
            self.vis = meshcat.Visualizer().open()
        self.vis.delete()


    def hooks(self):
        return {1:self.Hook1,
                2:self.Hook2,
                6:self.Hook6,
                7:self.Hook7}

    def dump(self):
        sleep(100)
        # wait forever

    @timeit
    @banner
    def Hook1(self, ice):
        "Draw the cell."
        logger = getLogger()
        v = self.vis["cell"]
        x,y,z = ice.repcell.mat
        count = 0
        for p,q,r in ((x,y,z),(y,z,x),(z,x,y)):
            for a in (np.zeros(3), p, q, p+q):
                draw_bond(v, f"cell{count}", a, r, 0.01, color=0x888888)
                count += 1


    @timeit
    @banner
    def Hook2(self, ice):
        "Output CoM of water molecules."
        logger = getLogger()
        if self.size_H > 0:
            return
        # prepare the reverse dict
        waters = defaultdict(dict)
        pos = ice.reppositions @ ice.repcell.mat
        v = self.vis["com"]
        for i, p in enumerate(pos):
            draw_atom(v, f"com{i}", p, 0.03, color=0x00ffff)
        return True


    @timeit
    @banner
    def Hook6(self, ice):
        "Output water molecules."
        logger = getLogger()
        logger.info("  Total number of atoms: {0}".format(len(ice.atoms)))
        # prepare the reverse dict
        waters = defaultdict(dict)
        molorder = -1
        for atom in ice.atoms:
            resno, resname, atomname, position, order = atom
            if resno == 0:
                molorder += 1
            if "O" in atomname:
                waters[molorder]["O"] = position
            elif "H" in atomname:
                if "H0" not in waters[molorder]:
                    waters[molorder]["H0"] = position
                else:
                    waters[molorder]["H1"] = position
        v = self.vis["water"]
        for i, water in waters.items():
            vv = v[f"water{i}"]
            O = water["O"]
            H0 = water["H0"]
            H1 = water["H1"]
            draw_atom(vv, f"O{i}", O, 0.03, color=0xff0000)
            draw_atom(vv, f"HA{i}", H0, self.size_H, color=0x00ffff)
            draw_atom(vv, f"HB{i}", H1, self.size_H, color=0x00ffff)
            draw_bond(vv, f"OHA{i}", O, H0-O, 0.01)
            draw_bond(vv, f"OHB{i}", O, H1-O, 0.01)
        v = self.vis["HB"]
        for i,j in ice.spacegraph.edges(data=False):
            if i in waters and j in waters:  # edge may connect to the dopant
                O = waters[j]["O"]
                H0 = waters[i]["H0"]
                H1 = waters[i]["H1"]
                d0 = H0 - O
                d1 = H1 - O
                rr0 = np.dot(d0,d0)
                rr1 = np.dot(d1,d1)
                if rr0 < rr1 and rr0 < 0.245**2:
                    draw_bond(v, f"HB{i}_{j}", O, d0, 0.005, color=0xffff00)
                if rr1 < rr0 and rr1 < 0.245**2:
                    draw_bond(v, f"HB{i}_{j}", O, d1, 0.005, color=0xffff00)
        self.nwateratoms = len(ice.atoms)


    @timeit
    @banner
    def Hook7(self, ice):
        "Output guest molecules. (WARNING: Not tested.)"

        def choose_palette(i):
            import colorsys
            w = (1+5**0.5)/2
            hue = (i/w) % 1
            r,g,b = colorsys.hsv_to_rgb(hue,1,1)
            return int(r*255)*65536+int(g*255)*255+int(b*255)

        logger = getLogger()
        logger.info("  Total number of atoms: {0}".format(len(ice.atoms)))
        gatoms = ice.atoms[self.nwateratoms:]
        palettes = dict()
        v = self.vis["guest"]
        # 原子列を、再度residueごとにグループ化する。無駄なように見えるが、
        # Stage7の内部処理がかなり複雑なので、Stage7結果をそのまま使うのが合理的。
        res_bucket = dict()
        molorder = -1
        for atom in gatoms:
            resno, resname, atomname, position, order = atom
            if resno == 0:
                molorder += 1
            reslabel = f"{resname}_{molorder}"
            if reslabel not in res_bucket:
                res_bucket[reslabel] = dict() # 原子名ごとの分類
            if atomname not in res_bucket[reslabel]:
                res_bucket[reslabel][atomname] = []
            res_bucket[reslabel][atomname].append(position)
        for reslabel in res_bucket:
            vv = v[reslabel]
            for atomname in res_bucket[reslabel]:
                if atomname not in palettes:
                    palettes[atomname] = choose_palette(len(palettes)+1)
                color = palettes[atomname]
                for i, position in enumerate(res_bucket[reslabel][atomname]):
                    draw_atom(vv, f"{atomname}_{i}", position, 0.04, color=color)
