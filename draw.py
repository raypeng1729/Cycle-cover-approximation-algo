import turtle

str1 = "test_data.txt"
str2 = "test_data_solution.txt"

with open(str1, "r") as inFile:
    n = int(inFile.readline())
    x = []
    y = []
    mx = 0
    for i in range(n):
        st = inFile.readline()
        num = st.split()
        num[0] = int(num[0])
        num[1] = int(num[1])
        x.append(num[0])
        y.append(num[1])
        mx = max(mx, num[0], -num[0], num[1], -num[1])

for i in range(n):
    x[i] = -x[i] / mx
    y[i] = -y[i] / mx

t = turtle.Turtle()
wn = turtle.Screen()
wn.setworldcoordinates(1.1, 1.1, -1.1, -1.1)
t.speed(10)
t.write("(0, 0)", font = ("Helvetica", 10))
t.goto(1, 0)
t.goto(-1, 0)
t.goto(0, 0)
t.goto(0, 1)
t.goto(0, -1)
t.up()
t.goto(1, 1)
ms = str(-mx)
t.write("(" + ms + ", " + ms + ")", font = ("Helvetica", 10))
t.down()
t.goto(1, -1)
t.goto(-1, -1)
s = str(mx)
t.write("(" + s + ", " + s + ")", font = ("Helvetica", 10))
t.goto(-1, 1)
t.goto(1, 1)
t.up()

ord = []
with open(str2, "r") as inFile:
    ord = inFile.readlines()
for i in range(n):
    ord[i] = int(ord[i]) - 1

t.goto(x[ord[0]], y[ord[0]])
t.down()

for i in range(n):
    t.dot(3)
    t.goto(x[ord[i]], y[ord[i]])
t.goto(x[ord[0]], y[ord[0]])

wn.exitonclick()
