import turtle as tt
from math import *
from time import sleep

# 定义颜色列表
colors = ["red", "blue", "green", "grey", "purple", "orange", "pink", "cyan"]

# 打开文件并读取对象个数
file = open("o9", "r")
num_objects = int(file.readline())

# 创建海龟对象列表
turtles = []
for i in range(num_objects):
    t = tt.Turtle("circle", 1)
    t.penup()
    # t.pensize(1)
    t.color(colors[i % len(colors)])  # 分配颜色
    t.resizemode("user")
    # t.shapesize(0.1, 0.1, 5)
    t.speed(0)
    turtles.append(t)

# 读取时间步数
time = int(file.readline())
for i in range(time):
    sleep(1/15)
    # print(file.readline(), end='')
    for t in turtles:
        x = float(file.readline())
        y = float(file.readline())
        t.goto(x, y)

tt.done()