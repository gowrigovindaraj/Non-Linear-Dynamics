# Visualisation of KAM Theory using Python
*****
Gowri Govindaraj

EP20BTECH11007

## Introduction
In the 1600s, Newton wrote down the set of differential equations that a system of massive bodies obey when they interact through the effects of gravitational force. It was observed that, if only 2 bodies were present, the system can be explicitly solved. It was also observed that the bodies moved along Keplerian ellipses around their collective Centre of Mass. 

For a more complex 3 body system, *no fixed, singular solution exists*. Consider our Solar System - as some bodies (planets), are much lighter than others (Sun). The mutual gravitational force between the planets is much weaker than any system where the Sun participates. 
Therefore, the system can be simplified:
- *Ignoring the interactions between the planets -> This gives an integrable system that can be solved as mentioned above*
- *Then we systematically include the interaction between the planets*

This method was extensively used throughout the 19<sup>th</sup> century, and series expansions were formed as a solution to these equations. However, as the terms had small denominators, their convergence could not be established.

If you have a typical nonlinear oscillator, you will get *resonances*  whenever the perturbing force has a frequency that is a rational multiple of the oscillator's natural frequency, since the nonlinearity will induce oscillations at all multiples of the fundamental driving frequency.
In a similar way, one planet exerts a periodic pull on the motion of another, and **if their orbital periods are equal, resonance and instability might occur**. In the perturbation theory, even if the two periods are not precisely equal, but merely roughly so, the effects cause *convergence issues*.

In 1954, **A.N.Kolmogorov** suggested a method, in the ICM, Amsterdam, to oversome these issues.
His suggestions contained two ideas which are central to all applications of the KAM techniques. These two basic ideas are:
- *Linearize the problem about an approximate solution and solve it.*
-  *Inductively improve the approximate solution by using the solution of the linearized problem as the basis of a Newton’s method argument.*

The proof of this was carried out by *Jürgen Moser*  and *Vladimir Arnold*, and hence this theory came to be known as KAM Theory (*Kolmogorov-Arnold-Moser theory*).
<div style="page-break-after: always;"></div>

## KAM Theory and its Mathematical Expression
Kolmogorov, in the conference, talked about the question of stability of the solar system. 
He proposed that despite the general chaotic behavior of the three–body problem, there could be “*islands of stability*” which were protected from chaos, allowing some orbits to remain regular even while other nearby orbits were highly chaotic.
While he came up with the brief outline on how to approach this, he could not complete the proof. 

Over the next 10 years, the work of German mathematician Jürgen Moser and Kolmogorov's former student Vladimir Arnold provided evidence of Kolmogorov's conjecture. The proof was based on a *series of integer ratios*  that approximated irrational numbers. KAM demonstrated that **by depending on the irrationality of the orbital period ratio, some orbits are shielded against adjacent chaos**.

#### Resonance in Ratios
Jupiter's orbital period is 11.86 years, however if it were exactly 12 years, it would have a 12:1 ratio with the Earth's orbital period. The term *resonance*  refers to the ratio of integers.
However, if this ratio were a low integer ratio, such as 4:3, the two planets would align every 12 years. This type of resonance with low integer ratios causes a large gravitational disturbance, which leads the orbit of the smaller planet to change. If the disturbance is large enough, it might cause the Earth's orbit to be disrupted, resulting in a chaotic course that could eventually dislodge the Earth from the solar system.

The planets have a hard time aligning as the resonance ratio reaches a ratio of large integers, such as 87:32, as KAM found, and the perturbation stays small. One intriguing aspect of this idea is that a close orbital ratio of 5:2 = 1.5, which is just slightly different from 87:32 = 1.7, might exist. The 5:2 resonance, on the other hand, may cause a lot of chaos, but the 87:32 resonance is essentially impervious to it. It is conceivable to have both chaotic and regular orbits in the same dynamical system in this fashion.

KAM Theory proposed the following:
> An irrational orbital ratio protects the regular orbits from chaos.

Therefore, we examine the most irrational number we know : The Golden Ratio.

#### Golden Ratio
**Golden ratio**, also known as the **golden section,** **golden mean**, or **divine proportion**, in mathematics,  denoted by the Greek letter ϕ or τ, is 
$$\phi = \frac{1+\sqrt{5}}{2} \approx 1.6180339887$$ 

We find the golden ratio when we divide a line into two parts so that: the whole length divided by the long part **_is also equal to_** the long part divided by the short part and is the hardest number to approximate with a ratio of small integers.

#### Solving the System using Python
As Kolmorogov proposed, because the dynamics of three-body systems are difficult to comprehend directly, *we assume that there are two large masses and a third small mass in a constrained three-body system*. The small mass has no effect on the dynamics of the two-body system in this way, therefore all we have to do is concentrate on the dynamics of the small body. 
This *reduces the dynamics to two dimensions* (the position and momentum of the third body), which makes visualisation much easier, but the dynamics still need differential equation solutions. The last step is to *replace the differential equations* with repeatedly solved simple difference equations.

The action-angle variables are connected by a perturbation in a simple discrete recurrent map that captures the key behaviour of the three-body issue. The Twist Map, the Chirikov Map, and the Standard Map are all variations on this form. The most important mapping is
> $$J_{n+1} = J_n + \epsilon sin{\theta_n}$$
$$\theta_{n+1} = mod(\theta_n + J_n, 2\pi)$$

Here, J represents the required Action variable, while $\theta$ represents the angle variable. 
$\epsilon$ is the perturbation parameter. 
The next step is to *initialise the above variables*, and the future values are *obtained by iteration*.
We run the python program for various values of $\epsilon$.

~~~py
import numpy as np

from scipy import integrate

from matplotlib import pyplot as plt

plt.close('all')

eps = 0.3 #setting a value of epsilon

np.random.seed(2) #generate most random numbers

plt.figure(1)

for iter in range(0,50):

    j_final = np.pi*(1.5*np.random.random()-0.5)  #action statement

    theta_final = 2*np.pi*np.random.random()         #angle statement

    orbit = np.int(200*(j_final+np.pi/2))  #initialising values

    rplot = np.zeros(shape=(orbit,))

    thetaplot = np.zeros(shape=(orbit,))

    x = np.zeros(shape=(orbit,))

    y = np.zeros(shape=(orbit,))    

    for iter2 in range(0,orbit):   #iterating

        j_new = j_final + eps*np.sin(theta_final)

        theta_new = np.mod(theta_final+j_new,2*np.pi)

        rplot[iter2] = j_new

        thetaplot[iter2] = np.mod(theta_new-np.pi,2*np.pi) - np.pi            

        j_final = j_new

        theta_final = theta_new

        x[iter2] = (j_new+np.pi+0.25)*np.cos(theta_new)

        y[iter2] = (j_new+np.pi+0.25)*np.sin(theta_new)

    plt.plot(x,y,'o',ms=1)  #making a plot for the various iterations

    plt.title("KAM Twist Map for ε = 0.3")

    plt.xlabel("J")

    plt.ylabel("θ")

plt.show()
~~~

If $\epsilon$ = 0; the perturbation = 0, and **all orbits are circular and regular**. 
![[eps = 0.png | 350]]

As the value of $\epsilon$ increases, **the open orbits breaks up into chains of closed orbits**. 
We observe the occurrence of **regular periodic, open orbits and regions of chaos**. 
![[eps = 0.3 1.png | 350]]

As the value of $\epsilon$ increases further,  so do the **regions exhibiting chaotic behaviour**( broken orbits).
![[eps = 0.8.png | 350]]

We observe the case with $\epsilon$ = 0.98, and see the regions where **regular orbits and chaotic orbits co-exist**. 
What makes such a condition is the orbital period ratio; for ratios that are sufficiently irrational, the regular orbits appear. For orbital ratios that consist of small integers, the perturbations caused drive the dynamic into chaos.


![[eps = 0.98 1.png | 350]]




## Application in our Solar System
Examine Kepler's Third Law in our solar system :
> The squares of the sidereal periods (of revolution) of the planets are directly proportional to the cubes of their mean distances from the Sun.

i.e., the square of the orbital period increases as the cube of the average radius

Consider the before mentioned 3 body system of the Earth, Jupiter and the Sun, where Earth is the smaller third body. 
Examining the stability of Earth's orbit as a function of its distance from the Sun, the orbital ratio wrt. Jupiter is smooth.

For Earth, this value is of the sorts 12:1 resonance, but as the distance from the Sun increases, the *ratio decreases*.
When the orbital ratio is adequately irrational, the orbit would be unaffected by Jupiter's gravitaional force. SImilarly, if the ratio is a smaller integral value, *perturbations occur and the effect increases*, and **regions of chaotic motion get vast**.
The regions of regular motion (orbital ratio irrational) act as a *barrier*, containing the range of chaotic orbits and restricting them.
In this way, number theory restrains the chaos of our solar system.

The asteroid belt provides a remarkable instance of the orbital resonance effect. The many tiny bodies operate as orbital resonance probes. 
The drag of Jupiter creates gaps in the distribution of asteroid radii, with large gaps, known as Kirkwood Gaps, appearing at orbital ratios of 3:1, 5:2, 7:3, and 2:1. The radii where chaotic activity occurs are these gaps, whereas the areas in between are stable. Because chaotic motion tends to push asteroids out of resonance zones, most asteroids spend the majority of their time in the stable regions. The same physics that causes gaps in Saturn's rings during resonances with its numerous moons is responsible for the Kirkwood gaps.




![[resonance-gap 1.jpg | 500]]

## Conclusion
The KAM theory may be extended to partial differential equations (PDEs) with a Hamiltonian structure in infinite dimension. The wave equation, the (stationary) Schrödinger equation, KdV, and others are examples of such equations. Nonlinear perturbations of these equations may be reduced to an unlimited number of connected dynamical (ordinary differential) equations under certain conditions (e.g., for the wave equation one obtains infinitely many coupled harmonic oscillators). The embedding of a linear quasi-periodic flow on a finite dimensions torus into the infinite dimensional phase space associated with the equation can then be found as quasi-periodic solutions. Motions that are almost-periodic have also been investigated (i.e., trajectories with infinitely many independent frequencies). Since the 1990s, several breakthroughs in these areas have been made.
*****
## References
- An Introduction to KAM Theory C. Eugene Wayne January 22, **2008**
- Arnold, V. I., From superpositions to KAM theory. _Vladimir Igorevich Arnold. Selected Papers_ **1997**
- A Lecture on the Classical KAM Theorem Jurgen Poschel
- The Kirkwood Gap image has been taken from [https://nineplanets.org/kirkwood-gap/](https://nineplanets.org/kirkwood-gap/)
- The KAM Story: A friendly introduction to the content, history and significance of Classical Kolmogorov-Arnold-Moser Theory.Dumas, H. S., World Scientific: **2014**.
- The python code used can be found at [<u>github hyperlink</u>](https://github.com/gowrigovindaraj/code/blob/main/ep20btech11007_nld.py)
*****