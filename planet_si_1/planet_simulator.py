import turtle as tt
from math import *
from time import sleep

p0 = tt.Turtle("circle", 1)
p1 = tt.Turtle("circle", 1)
p2 = tt.Turtle("circle", 1)
p3 = tt.Turtle("circle", 1)
p0.penup()
#p1.penup()
p2.penup()
p3.penup()
p0.speed(0)
p1.speed(0)
p2.speed(0)
p3.speed(0)

# with open("/Volumes/Untitled/cppworkspace/planet_simulator/orbit", "r") as file:
file =  open("./6", "r")
time = int(file.readline())
for i in range(time):
    sleep(0.001)
    gaps = file.readline()
#    print(gaps, end='')
    p0x = float(file.readline()) / 1e9
    p0y = float(file.readline()) / 1e9
    p1x = float(file.readline()) / 1e9
    p1y = float(file.readline()) / 1e9
    p2x = float(file.readline()) / 1e9
    p2y = float(file.readline()) / 1e9
    p3x = float(file.readline()) / 1e9
    p3y = float(file.readline()) / 1e9
    p0.goto(p0x, p0y)
    p1.goto(p1x, p1y)
    p2.goto(p2x, p2y)
    p3.goto(p3x, p3y)

tt.done()
