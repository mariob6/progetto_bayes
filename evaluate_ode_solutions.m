y0 =[7;-10];
tspan = [0:0.5:60];
f1 = @(t,y) [72./(36 + y(2))-1; y(1)-1];

[T,Y] = ode45(f1,tspan,y0);

csvwrite('toy_example_sol.csv',Y)