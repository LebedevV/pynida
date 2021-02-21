#!/usr/bin/env python

import sys,os
from PyQt5 import QtWidgets
from PyQt5.QtGui import QFont

from pynida.launch import mywindow

import pynida.classes

def main():
	app = QtWidgets.QApplication([])

	app.setFont(QFont('Arial',10))

	application = mywindow()
	application.show()
	sys.exit(app.exec())
	
if __name__ == "__main__":
	main
