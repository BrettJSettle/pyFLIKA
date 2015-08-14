echo off
IF "%1"=="-install" (
	call pip install -r list.txt > install.txt
	del install.txt
)
python pyFlika.py -no3D
del /S *.pyc
