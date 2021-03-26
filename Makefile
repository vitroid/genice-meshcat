.DELETE_ON_ERROR:
GENICE=genice2
BASE=genice2_meshcat
PACKAGE=genice2-meshcat

all: README.md

%: temp_% replacer.py $(wildcard $(BASE)/formats/*.py) $(wildcard $(BASE)/*.py)
	pip install genice2_dev svgwrite
	python replacer.py < $< > $@

test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install pillow
	pip install --index-url https://test.pypi.org/simple/ $(PACKAGE)

install:
	python3 setup.py install
uninstall:
	-pip uninstall -y $(PACKAGE)
build: README.md $(wildcard $(BASE)/lattices/*.py) $(wildcard $(BASE)/*.py)
	python ./setup.py sdist bdist_wheel

deploy: build
	twine upload dist/*
check:
	./setup.py check
clean:
	-rm $(ALL) *~ */*~
	-rm -rf build dist *.egg-info
	-find . -name __pycache__ | xargs rm -rf
