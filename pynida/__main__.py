#!/usr/bin/env python

import sys
import os
from PyQt5 import QtWidgets
from PyQt5.QtGui import QFont

sys.path.append(os.path.abspath('.'))

from pynida.launch import mywindow
from pynida.config import GlobalFont


def main():
	app = QtWidgets.QApplication([])

	app.setFont(QFont(GlobalFont.NAME, GlobalFont.SIZE))

	application = mywindow()
	application.show()
	sys.exit(app.exec())


if __name__ == "__main__":
	main()
