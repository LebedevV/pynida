# pynida

This software allows one to perform a full cycle of *in-situ* nanoindentation data analysis, including the video processing

**Current limitations:**
---
- Only *avi* files for video  
- Dark-field images are expected for the particle shape evolution analysis
- Magnification should not be shanged during the experiment

**Installation**
---

### Requirements
matplotlib  
numpy  
opencv-python  
pandas  
PyQt5  
scipy  
\>=scikit-image-0.18

##### Currently, to install the software one should install all the prerequisites, and to execute in the folder for choise the following commands:
*\$git init*  
*\$git pull https://github.com/LebedevV/pynida*
##### Afterwards, to run it:
*\$python pynida*

### Feature requests and planned work

**Known bugs**
- fonts in GUI might be corrupted 
- re-measuring of the rotation angle happens during the comparison stage. Might results in the wrong rotation angle in output
- (not related to the software itself) scales measured by the scalebar and by metadata may be different
-	(to check) button "q" for exit from plt windows is not working on Ubuntu
- (to check) test imgs are saving nearby the opened file, not in the output dir
 
 **Planned work, long-term**
- add dm3/dm3 import
- scale down preview (need in a lineEdit for that?)
- add "save video" option
- fully automatic scale detection. Will allow one to deal with the variable magnification
- for rotation: GUI with grid and drag&drop
- (?) Graphical scheme of indentation to GUI
 
 **Planned work, short-term**
- proper logs should be configured  
- allow shape analysis for the BF-TEM  
- (partially done) preview of areas/contours, in pix/fr
- Fix column names in output files  
- (partially done) Add descriptions to output files
- (partially done) file format checks
- update *pnm* by units recall
- make code more readable, especially in the areas analysis
- (done) probe the write access to the output dir
- reduce number of global variables
- (almost done) move errors to qmesg
- same checks on different buttons - to unify
- video is probing every time in update_crop. To replace with the affine transform
- create a list of failed frames for areas/contours
- add unit tests!
  
