#!/usr/bin/env python

import sys,os
from configparser import ConfigParser
from PyQt5 import QtWidgets
from PyQt5.QtGui import QFont

from pynida.launch import mywindow

import pynida.classes

def main():
	app = QtWidgets.QApplication([])

	config = ConfigParser()
	config.read('pynida/config.cfg')

	app.setFont(
		QFont(config["DEFAULT"]["GLOBAL_FONT_NAME"],
			  int(config["DEFAULT"]["GLOBAL_FONT_SIZE"])))

	application = mywindow()
	application.show()
	sys.exit(app.exec())
	
if __name__ == "__main__":
	main
