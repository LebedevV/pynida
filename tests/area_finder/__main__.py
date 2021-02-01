print("____-menu_____") 
print("1: to find area of square \n\ 
2: to find area of rectangle\n\ 
3: to find area of circle") 

ch = int(input()) 

if ch == 1: 
	from pynida.square import square
	print("enter side") 
	s = int(input()) 
	print ("the area is ", square(s)) 

if ch == 2: 
	from pynida.rectangle import rectangle
	print("enter length and breadth") 
	l = int(input()) 
	b = int(input()) 
	print("the area is ", rectangle(l, b)) 

if ch == 3: 
	from pynida.circle import circle
	print("enter radius") 
	r = int(input()) 
	print("the area is ", circle(r)) 

