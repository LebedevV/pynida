#!/usr/bin/env python

import sys
from configparser import ConfigParser
from PyQt5 import QtWidgets
from pynida.launch import mywindow

if __name__ == "__main__":
    app = QtWidgets.QApplication([])

    config = ConfigParser()
    config.read('pynida/config.cfg')

    application = mywindow(config)
    application.show()
    sys.exit(app.exec())
