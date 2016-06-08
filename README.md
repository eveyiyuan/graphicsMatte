This project implements a automatic 2D image matte

This project relies on the OpenCV library and the Figtree Library
(vmorariu/figtree)

In the include folder, please include Figtree's figtree.h,
figtree_internal.h, and KCenterClustering.h

In the src folder, please include Figtree's figtree.cpp and
KCenterClustering.cpp

# Instructions:
* to make: make
* to run: ./matting /<path to image to matte/> /<path to new background image/>
* interactive:
	1. When your first image appears, hold down the right mouse button to draw white scribbles to denote the foreground
	2. Then, hold down the left mouse button to draw black scribbles to denote the background
	3. Then, press space to begin the matte calculation
	4. A Black and White image of the matte will then appear
	5. Press space to begin the Poisson image editing portion
	6. Your second image will appear
	7. Right click on the image to indicate where you want the foreground from the previous image to be centered
	8. Hold down the left mouse button and draw a circle around where you right-clicked to indicate the size you want the foreground to be 
	9. Press space to generate your merged image
