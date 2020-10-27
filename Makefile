install:
	git clone https://github.com/RussTedrake/meshcat-python --branch colab --recursive
	cd meshcat-python; python3 setup.py install
	python3 setup.py install
