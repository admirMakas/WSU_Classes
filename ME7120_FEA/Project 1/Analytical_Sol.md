

```python
import sympy as sp
sp.init_printing(use_latex='mathjax')
from IPython.display import display
import numpy as np

kn = sp.Symbol('k_n', real = True)
lamb = sp.Symbol('lambda', real = True)
omegan = sp.Symbol('omega_n', real = True)
E = sp.Symbol('E', real = True)
I = sp.Symbol('I', real = True)
L = sp.Symbol('L', real = True)
w = sp.Symbol('w', real = True)
x = sp.Symbol('x', real = True)
q = sp.Symbol('q', real = True)
t = sp.Symbol('t', real = True)
X = sp.Symbol('X', real = True)
f = sp.Symbol('f', real = True)
A1 = sp.Symbol('A_1', real = True)
A2 = sp.Symbol('A_2', real = True)
A3 = sp.Symbol('A_3', real = True)
A4 = sp.Symbol('A_4', real = True)
```

For a homogeneous Euler-Bernoulli beam the equation of motion is:


```python
E*I*sp.diff(w(x, t), x, 4) + lamb*sp.diff(w(x, t), t, 2) - q
```




$$E I \frac{\partial^{4}}{\partial x^{4}}  w{\left (x,t \right )} + \lambda \frac{\partial^{2}}{\partial t^{2}}  w{\left (x,t \right )} - q$$



The free vibration case reduces to:


```python
Eq = E*I*sp.diff(w(x, t), x, 4) + lamb*sp.diff(w(x, t), t, 2)
Eq
```




$$E I \frac{\partial^{4}}{\partial x^{4}}  w{\left (x,t \right )} + \lambda \frac{\partial^{2}}{\partial t^{2}}  w{\left (x,t \right )}$$



Assume following form of $w(x,t) = X(x)f(t)$


```python
Eq = E*I*sp.diff(X(x)*f(t), x, 4) + lamb*sp.diff(X(x)*f(t), t, 2)
Eq
```




$$E I f{\left (t \right )} \frac{d^{4}}{d x^{4}}  X{\left (x \right )} + \lambda X{\left (x \right )} \frac{d^{2}}{d t^{2}}  f{\left (t \right )}$$



Divide above expression by $\lambda X(x) f(t)$


```python
Eq1 = sp.nsimplify(Eq/(lamb*X(x)*f(t)))
Eq1.expand()
```




$$\frac{E I \frac{d^{4}}{d x^{4}}  X{\left (x \right )}}{\lambda X{\left (x \right )}} + \frac{\frac{d^{2}}{d t^{2}}  f{\left (t \right )}}{f{\left (t \right )}}$$



By seperation of variables each of the terms is constant with respect to one another. This constant is defined as $\omega_n^2$, which is the natural frequency.


```python
Eqx = (Eq1.expand().coeff((E*I)/lamb)*((E*I)/lamb)-omegan**2)*((lamb*X(x))/(E*I))
Eqx.expand()
```




$$\frac{d^{4}}{d x^{4}}  X{\left (x \right )} - \frac{\lambda \omega_{n}^{2}}{E I} X{\left (x \right )}$$



Above experssion in terms of $x$ can be used to solve for natural frequencies and mode shapes. Following expression satisfies the above PDE.


```python
w = A1*sp.sin(kn*x) + A2*sp.cos(kn*x) + A3*sp.sinh(kn*x) + A4*sp.cosh(kn*x)
w
```




$$A_{1} \sin{\left (k_{n} x \right )} + A_{2} \cos{\left (k_{n} x \right )} + A_{3} \sinh{\left (k_{n} x \right )} + A_{4} \cosh{\left (k_{n} x \right )}$$



Where $k_n^4 = \frac{\lambda \omega_n^2}{EI}$

For cantilevered beam, following B.C.'s are used:
- $w\ =\ 0$ at $x=0$, zero displacement  
- $\frac{dw}{dx} = 0$ at $x=0$, zero slope
- $\frac{d^2w}{dx^2} = 0$ at $x=L$, zero bending moment
- $\frac{d^3w}{dx^3} = 0$ at $x=L$, zero shear

For first B.C. $w\ =\ 0$ at $x=0$


```python
w1 = w.subs(x, 0)
w1
```




$$A_{2} + A_{4}$$



For second B.C. $\frac{dw}{dx} = 0$ at $x=0$


```python
w2 = sp.diff(w, x).subs(x, 0)
w2
```




$$A_{1} k_{n} + A_{3} k_{n}$$



For third B.C. $\frac{d^2w}{dx^2} = 0$ at $x=L$


```python
w3 = sp.diff(w, x, 2).subs(x, L)
w3.expand()
```




$$- A_{1} k_{n}^{2} \sin{\left (L k_{n} \right )} - A_{2} k_{n}^{2} \cos{\left (L k_{n} \right )} + A_{3} k_{n}^{2} \sinh{\left (L k_{n} \right )} + A_{4} k_{n}^{2} \cosh{\left (L k_{n} \right )}$$



Sub in following relationships and solve for $A_3$:
- $A_1$ = $-A_3$
- $A_2$ = $-A_4$


```python
w3s = sp.solve(w3.expand().subs(A1, -A3).subs(A2, -A4), A3)
A3s=w3s[0]
A3s
```




$$- \frac{A_{4} \left(\cos{\left (L k_{n} \right )} + \cosh{\left (L k_{n} \right )}\right)}{\sin{\left (L k_{n} \right )} + \sinh{\left (L k_{n} \right )}}$$



For fourth B.C. $\frac{d^3w}{dx^3} = 0$ at $x=L$


```python
w4 = sp.diff(w, x, 3).subs(x, L)
w4.expand()
```




$$- A_{1} k_{n}^{3} \cos{\left (L k_{n} \right )} + A_{2} k_{n}^{3} \sin{\left (L k_{n} \right )} + A_{3} k_{n}^{3} \cosh{\left (L k_{n} \right )} + A_{4} k_{n}^{3} \sinh{\left (L k_{n} \right )}$$



Sub in following relationships and solve for $A_4$:
- $A_1$ = $-A_3$
- $A_2$ = $-A_4$


```python
w4s = sp.solve(w4.expand().subs(A1, -A3).subs(A2, -A4), A4)
A4s=w4s[0]
A4s
```




$$\frac{A_{3} \left(\cos{\left (L k_{n} \right )} + \cosh{\left (L k_{n} \right )}\right)}{\sin{\left (L k_{n} \right )} - \sinh{\left (L k_{n} \right )}}$$



**In next step, expression derived above for $A_1,\ A_2$ and $A_3$ can be plugged into equation for $A_4$ ($4^{th}$ boundary condition) to get following term.**


```python
test = w4.subs(A3, A3s).subs(A1, -A3s).subs(A2, -A4)
test.simplify()
```




$$- \frac{2 A_{4} k_{n}^{3} \left(\cos{\left (L k_{n} \right )} \cosh{\left (L k_{n} \right )} + 1\right)}{\sin{\left (L k_{n} \right )} + \sinh{\left (L k_{n} \right )}}$$



**From this term it is evident that non trivial solutions exist only if $cos(Lk_n)cosh(Lk_n) + 1 = 0$. This equation can be solved numerically for roots of $Lk_n$, which will yield the natural frequencies.**  

**First 10 solutions are numerically determined below:**


```python
x2=np.linspace(0, 30, 100000)

#print(x)

Y=np.cos(x2)*np.cosh(x2)+1

S = len(Y)
knL = []
for i in range(S-1):
    if Y[i]<0.0 and Y[i-1]>0:
        #print('First Condition')
        #print('Fun Values =', Y[i], Y[i-1])
        #print('Wave Num Values =', x[i], x[i-1])
        value = (x2[i] + x2[i-1])/2
        #print(value)
        knL.append(value)
    elif Y[i]>0 and Y[i-1]<0:
        #print('Second Condition')
        #print('Fun Values =', Y[i], Y[i-1])
        #print('Wave Num Values =', x[i], x[i-1])
        value = (x2[i] + x2[i-1])/2
        #print(value)
        knL.append(value)
print(knL)
```

    [1.8751687516875166, 4.6939969399693995, 7.8548285482854823, 10.995559955599555, 14.137191371913719, 17.278822788227881, 20.420454204542043, 23.562085620856209, 26.703417034170339, 29.845048450484505]
    

**For a cantilever beam with following properties, the first 10 natural frequencies in Hz are calculated below using following expression $\omega_n = k_n^2 \sqrt{\frac{EI}{\rho A}}$, since $k_n L = x$, the formula can be expressed as $\omega_n = \frac{x^2}{L^2} \sqrt{\frac{EI}{\rho A}}$**  

- L = 25
- b = 0.02
- h = 0.02
- $\rho$ = 7850 $\frac{kg}{m^3}$ (Steel Density)
- E = $2.03^{11}$ Pa


```python
L = 25
b = 0.02
h = 0.02
rho = 7850
E = 2.03e11
A = b*h
Ix = (b*h**3)/12

for j in range(len(knL)):
    omega = ((knL[j]**2/L**2)*sp.sqrt((E*Ix)/(rho*A)))/(2*np.pi)

    print(omega)
```

    0.0262889591712005
    0.164732117817639
    0.461281563745762
    0.903914554844224
    1.49423516054184
    2.23213767928907
    3.11762211108591
    4.15068845593235
    5.33121692466034
    6.65943300260438
    

### **Next step is to obtain mode shapes for the above requencies.**

**Plugging in expressions for $A_1,\ A_2,\ A_3,\ A_4$ into the assumed solution for $w$ based on the prescribed boundary conditions yields the equation for the mode shapes**


```python
ws1 = (-A3s*sp.sin(kn*x)+A3s*sp.sinh(kn*x)).factor()
ws2 = -A4*sp.cos(kn*x)+A4*sp.cosh(kn*x)
#display(ws1)
#display(ws2)
ws = (ws1+ws2).simplify()
display(ws)
```


$$- \frac{A_{4}}{\sin{\left (L k_{n} \right )} + \sinh{\left (L k_{n} \right )}} \left(\left(\sin{\left (L k_{n} \right )} + \sinh{\left (L k_{n} \right )}\right) \left(\cos{\left (k_{n} x \right )} - \cosh{\left (k_{n} x \right )}\right) - \left(\sin{\left (k_{n} x \right )} - \sinh{\left (k_{n} x \right )}\right) \left(\cos{\left (L k_{n} \right )} + \cosh{\left (L k_{n} \right )}\right)\right)$$


**Constant $A_4$, which is unique for each frequency, is considered arbitrary in this context since we are only after the mode shape. To plot mode shapes it is customary to use $A_4$ = $\frac{1}{2}$.**


```python
import matplotlib.pyplot as plt
import numpy as np
#define beam lenght as array
x1=np.linspace(0, 25, 5000)

# constant for defining mode
m = 0


# calculate k1
k1 = knL[m]/L

# define constant
c=0.5

# equation for mode shapes
mode1 = (c/(np.sin(knL[m]) + np.sinh(knL[m])))*\
        ((np.sin(k1*x1) - np.sinh(k1*x1))*(np.cos(knL[m]) + np.cosh(knL[m])))\
        -c*np.cos(k1*x1) + c*np.cosh(k1*x1)

#Plot of numerical and analytical solutions
fig = plt.figure(figsize = (15,8))

plt.plot(x1, mode1, 'b-', label = 'Mode = 1', linewidth = 4) #FFT solution plot

plt.legend(loc = 'upper right')
fig.suptitle('Mode Shapes', fontsize = 14)
plt.xlabel('Distance from end')
plt.ylabel('Displacement $w$, (arbitrary)')
plt.axis([-0.2, 25.2, -1.5, 1.5])

#plt.grid()
plt.show()
```


```python
from IPython.display import HTML

HTML('''<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
The raw code for this IPython notebook is by default hidden for easier reading.
To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>.''')
```




<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
The raw code for this IPython notebook is by default hidden for easier reading.
To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>.


