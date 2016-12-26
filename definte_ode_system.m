function dydt = odefun(t,y)
dydt = [72/(36 + y(2)) -1; y(1) -1]

