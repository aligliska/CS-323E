# Alice Gee, Mohammad Aga, Andrew Yang
# eid: ag67642, mba929, ay6764

def taylor_series(x):
    status = True
    multiplier = 1
    n = 1
    store = 0
    while status == True:
        temp = multiplier * (1/n) * (x**n)
        store += temp 
        n += 2
        multiplier *= -1

        if round(4*store, 6) == 3.141592:
            status = False
            
    return(round(4*store, 6), round(n/2))

print("The target value was:", taylor_series(1)[0])
print(f"It took approximately {taylor_series(1)[1]} loops to reach this value.")
print()


multiplier = 1
x = 1
num = 900000
store = 0
for i in range(1, num*2, 2):
    temp = multiplier * (1/i) * (x**i)
    store += temp
    multiplier *= -1
print("Manual version to get value:", round(4*(store), 6))
