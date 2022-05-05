if False:
    from scipy.integrate import odeint
    y0 = [50.0, 0.0,0.0]
    z = 10
    K = 1
    b = 5
    time_l = 10
    def func(y, t, z, K,b):
        xy,x,y = y
        dydt = [-K*xy*(z/(z+b)),K*xy*(z/(z+b)), K*xy*(z/(z+b))]
        return dydt
    t = np.linspace(0, time_l, 101)
    sol = odeint(func, y0, t, args=(z,K,b))
    xy = [sol[i][0] for i in range(len(sol))]
    x = [sol[i][1] for i in range(len(sol))]
    y = [sol[i][2] for i in range(len(sol))]
    plot(x,'x')