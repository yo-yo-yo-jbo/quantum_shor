# Quantum computing: QFT and Shor's algorithm
In this blogpost I'd like to describe how Shor's algorithm works, and how we're able to efficiently factorize large numbers on a Quantum computer.   
We will cover several mathematical topics that come from number theory, exponents of complex numbers and Fourier transforms.  
I highly advise to read my previous blogposts on the topic: [Introduction to Quantum computing](https://github.com/yo-yo-yo-jbo/quantum_computing_intro/) and [Grover's algorithm](https://github.com/yo-yo-yo-jbo/quantum_computing_grover/).

## About factorization and non-trivial divisors
The problem of [factorization](https://en.wikipedia.org/wiki/Factorization) is something you learn as a child: given a number, represent it as the product of prime numbers.  
Of course, small numbers can be factored "by hand", e.g. $1960 = 2^3 \cdot 5 \cdot 7^2$.  
The first thing you might notice is that finding a prime divisor of the target number reduces the problem. In our case, if you notice 1960 is even then you divide it by 2 (division is a relatively trivial operation) and get 980, and we have reduced the problem.  
The second thing (less trivial maybe) is that finding numbers that are not co-prime with the target number is also helpful. If we find a number $x$ with $1 \lt gcd \left( x, 1960 \right) \lt 1960$, then $d = gcd \left( x, 1960 \right)$ splits our number into $d$ and $\frac{1960}{d}$, both of which are strictly smaller than 1960, and so, factoring them becomes an easier task, and can be done recursively.  
Therefore, finding any non-trivial divisor of a number leads to solving the problem of factorization.

## The relation to function period finding
Larger numbers require computers and might use much more sophisticated algorithms (e.g. [Pollard's Rho algorithm](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm), [Quadratic sieve](https://en.wikipedia.org/wiki/Quadratic_sieve)) but none of them are practical when it comes to large numbers, such as the ones generated in RSA (which I [wrote about](https://github.com/yo-yo-yo-jbo/rsa_math/)).  
The Quadratic Sieve, for instance, attempts to find numbers $x$ and $y$ that yield $x^2 \equiv y^2 \pmod{n}$ but $x \not\equiv \pm y \pmod{n}$. Note that $x^2 - y^2 = (x-y)(x+y)$, and so $n$ divides the product $(x-y)(x+y)$ while dividing neither factor on its own, which is exactly what $x \not\equiv \pm y \pmod{n}$ tells us. That means that $gcd \left( x-y, n \right)$ is a non-trivial divisor of $n$.  
What Shor's algorithm will do is similar, but with $y=1$. To do that, we can try to find the period of the function $a^x \pmod{n}$; let's unpack what that means.  
Before we start, we assume $n$ is odd and is not a prime power. These cases are easy to detect and to handle classically: if $n$ is even we simply divide by 2, if $n$ is prime a [primality test](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test) tells us so and we are done, and if $n = p^k$ for $k \gt 1$ we can find $p$ by taking integer $k$-th roots for the few relevant values of $k$. This matters, since for those cases the argument below cannot produce a split:
1. Let us assume some number $1 \lt a \lt n - 1$. We also assume $gcd \left( a, n \right) = 1$, since otherwise $gcd \left( a, n \right)$ is already a non-trivial factor of $n$ and we got lucky.
2. If we multiply $a$ by itself over and over again we will eventually get to $1 \pmod{n}$. The reason for that is that there are only finitely many values modulo $n$, so the powers of $a$ must eventually repeat, and since $a$ is invertible modulo $n$ we can cancel and land back on 1. (Concretely, [Euler's theorem](https://en.wikipedia.org/wiki/Euler%27s_theorem) tells us $a^{\varphi \left( n \right)} \equiv 1 \pmod{n}$). We call the minimal number $r \gt 0$ which yields $a^r \equiv 1 \pmod{n}$ **the order** of $a \pmod{n}$. Let us assume we magically know that value.
3. Let us assume that order is even. Then: $a^r - 1 = \left( a^{\frac{r}{2}} \right)^2 - 1 = \left( a^{\frac{r}{2}} - 1 \right) \left( a^{\frac{r}{2}} + 1 \right) \equiv 0 \pmod{n}$.
4. Let us note $a^{\frac{r}{2}} - 1 \not\equiv 0 \pmod{n}$. If it were, then $\frac{r}{2}$ would be positive but smaller than the order $r$, which is a contradiction to the fact the order is the minimal positive value that yields $1 \pmod{n}$.
5. If $a^{\frac{r}{2}} + 1 \not\equiv 0 \pmod{n}$, or, in other words, $a^{\frac{r}{2}} \not\equiv -1 \pmod{n}$, then we have two numbers multiplied together that are non-zero and yield $0 \pmod{n}$. That means we get non-trivial divisors from them by $gcd \left( a^{\frac{r}{2}} + 1, n \right)$ and $gcd \left( a^{\frac{r}{2}} - 1, n \right)$.

That is the entire idea behind the classical part of Shor's algorithm, and it depends on one algorithmic assumption: find the order of $a \pmod{n}$.  
Note that I had other assumptions, such as the order being even and $a^{\frac{r}{2}} \not\equiv -1 \pmod{n}$. For an odd $n$ that is not a prime power, a randomly chosen $a$ satisfies both conditions with probability of at least $\frac{1}{2}$, so a handful of attempts suffices. If an $a$ turns out to be bad, we simply choose a different candidate $a$ randomly and try again. We will see later that the Quantum part can also fail occasionally, and that there are cheap ways of recovering from that without starting from scratch.

Note in the title I wrote "function period finding". A **periodic function** is a function that repeats itself after some time. Formally, we say that a function has a period of $p$ if $f \left( x \right) = f \left( x + p \right)$ for every $x$, and we are normally interested in the *minimal* such positive $p$. In our case for Shor's algorithm:
- The function is $f \left( x \right) = a^x \pmod{n}$.
- The period that we are looking for is exactly the order $r$ of $a \pmod{n}$. Indeed, $a^{x+r} = a^x \cdot a^r \equiv a^x \pmod{n}$.

So, if we had a magical algorithm to assess the period of a given function efficiently, we would know how to factor large numbers with relative ease.

## Background: DFT
Let us take a break from number theory and talk about a completely different topic: the [Discrete Fourier Transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform) (DFT).  
The general idea of a Fourier transform is to take a signal and describe it as a sum of "pure frequencies". Instead of asking "what is the value of the signal at time $x$" we ask "how much of frequency $k$ does the signal contain". Since periodicity is exactly a statement about frequency, this is a very natural tool for us; we are, after all, hunting for the period of $a^x \pmod{n}$.

The discrete version works on a finite vector of $N$ numbers $\left( x_0, x_1, \ldots, x_{N-1} \right)$ and produces a new vector of $N$ numbers $\left( y_0, y_1, \ldots, y_{N-1} \right)$ defined as:

$$y_k = \frac{1}{\sqrt{N}} \sum_{j=0}^{N-1} x_j \cdot e^{\frac{2 \pi i j k}{N}}$$

(Some textbooks write $e^{-\frac{2 \pi i j k}{N}}$ here. The two conventions differ only by which of the transform and its inverse you call "the" DFT. I am using the sign that matches the Quantum version we are heading towards.)  
The strange looking $e^{\frac{2 \pi i j k}{N}}$ is a [root of unity](https://en.wikipedia.org/wiki/Root_of_unity): if we denote $\omega = e^{\frac{2 \pi i}{N}}$ then $\omega^N = 1$, and the numbers $1, \omega, \omega^2, \ldots, \omega^{N-1}$ are $N$ points spread evenly on the unit circle in the complex plane. The $\frac{1}{\sqrt{N}}$ in front is just a normalization factor. It makes the transform preserve lengths, which will turn out to be important for us.

## Why the DFT finds periods
The entire magic of the DFT comes from one small fact about roots of unity: summing all of them gives zero. Indeed, for any integer $k$ that is not a multiple of $N$, the sum $\sum_{j=0}^{N-1} \omega^{jk}$ is a [geometric series](https://en.wikipedia.org/wiki/Geometric_series) with ratio $\omega^k \neq 1$, and therefore equals $\frac{\omega^{Nk} - 1}{\omega^k - 1} = \frac{1 - 1}{\omega^k - 1} = 0$. When $k$ *is* a multiple of $N$ every term is simply 1 and the sum is $N$. Intuitively, evenly spread vectors on the circle cancel each other out, unless they all point in the same direction.

This is exactly the "detector" we want:
- Suppose our input vector is periodic with period $r$, say $x_j = 1$ whenever $j$ is a multiple of $r$, and $x_j = 0$ otherwise.
- For a frequency $k$ that "matches" the period (that is, $\frac{kr}{N}$ is an integer) all the terms $\omega^{jk}$ point in the same direction and add up constructively, so $y_k$ becomes large.
- For any other frequency the terms are spread around the circle and cancel out, so $y_k$ becomes tiny.

In other words, the DFT takes something that is periodic in the *value* domain and turns it into something concentrated in the *frequency* domain: a needle instead of a haystack. Note the pleasant duality here: a signal with period $r$ produces spikes spaced $\frac{N}{r}$ apart, so reading the spacing of the spikes tells us $r$.  
There is one subtlety I swept under the rug: this picture is exact only when $r$ divides $N$. In general $\frac{N}{r}$ is not an integer, so the spikes sit *between* integer frequencies, and what we get are peaks that are slightly smeared over the nearest integers. They are still very much peaks, and we will see later how to deal with the smearing.

It is also worth writing the DFT as a matrix: the matrix $F$ whose entry in row $k$ and column $j$ is $\frac{1}{\sqrt{N}} \omega^{jk}$. With the normalization we chose, $F$ is a [unitary matrix](https://en.wikipedia.org/wiki/Unitary_matrix), i.e. $F F^{\dagger} = I$, which means it is a reversible, length-preserving linear operation. Remember that, it is about to become the whole point.

The catch is efficiency. Computing all $N$ outputs directly costs $N$ multiplications each, so $O \left( N^2 \right)$ operations. The celebrated [Fast Fourier Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform) (FFT) brings this down to $O \left( N \log N \right)$, which is excellent. However, for our purposes $N$ will be around $n^2$ (we will see exactly why later), i.e. exponential in the number of *bits* of $n$. Even writing down the input vector is hopeless. This is where the Quantum computer walks in.

## Background: QFT
Here is the observation that starts everything: a Quantum state of $q$ qubits is *already* a vector of $N = 2^q$ complex numbers (the amplitudes). And a Quantum gate is *exactly* a unitary matrix acting on that vector. We just said the DFT matrix is unitary, so the DFT is a legal Quantum gate. That gate is called the **Quantum Fourier Transform** (QFT).  
(Note I use $q$ for the number of qubits, since $n$ is already taken by the number we wish to factor.)

Formally, the QFT acts on basis states as:

$$QFT \left( \left| j \right\rangle \right) = \frac{1}{\sqrt{N}} \sum_{k=0}^{N-1} e^{\frac{2 \pi i j k}{N}} \left| k \right\rangle$$

And by linearity, on a general state $\sum_j x_j \left| j \right\rangle$ it produces $\sum_k y_k \left| k \right\rangle$ where the $y_k$ are precisely the DFT coefficients from before. Nothing about the mathematics changed, only *where the numbers live*. Instead of $N$ numbers sitting in memory, they are $N$ amplitudes riding on $q = \log_2 N$ qubits.

### Building the circuit
It is not obvious that such a huge matrix can be applied with few gates, so let us see why it can. Write $j$ in binary as $j = j_1 j_2 \ldots j_q$, meaning $j = \sum_{l=1}^{q} j_l 2^{q-l}$, and use the standard binary fraction notation $0.j_l j_{l+1} \ldots j_q = \frac{j_l}{2} + \frac{j_{l+1}}{4} + \ldots + \frac{j_q}{2^{q-l+1}}$.  
With a bit of algebra, the sum above factors into a product, the so-called *product representation*:

$$QFT \left( \left| j \right\rangle \right) = \frac{1}{\sqrt{N}} \bigotimes_{l=1}^{q} \left( \left| 0 \right\rangle + e^{2 \pi i \cdot 0.j_{q-l+1} \ldots j_q} \left| 1 \right\rangle \right)$$

Why does it factor? Because $e^{\frac{2 \pi i j k}{N}}$ splits over the binary digits of $k$, and once summed, each qubit of the output ends up depending only on a phase; it is never entangled with the others. Also note that $e^{2 \pi i m} = 1$ for any integer $m$, which is why only the *fractional* part of $\frac{j}{2^l}$ survives in each factor.

This gives us a recipe, since each output qubit is just $\left| 0 \right\rangle$ plus $\left| 1 \right\rangle$ with some rotation:
1. Apply a [Hadamard gate](https://en.wikipedia.org/wiki/Quantum_logic_gate#Hadamard_gate) to the first qubit. This produces $\left| 0 \right\rangle + e^{2 \pi i \cdot 0.j_1} \left| 1 \right\rangle$, since the Hadamard adds a minus sign (which is $e^{\pi i}$) exactly when $j_1 = 1$.
2. Add the remaining bits of the phase with controlled rotations. The gate $R_m$ leaves $\left| 0 \right\rangle$ untouched and multiplies $\left| 1 \right\rangle$ by $e^{\frac{2 \pi i}{2^m}}$. Applying $R_2$ controlled by qubit 2, then $R_3$ controlled by qubit 3, and so on up to $R_q$, appends the digits $j_2, j_3, \ldots, j_q$ to the binary fraction one at a time. A controlled rotation only fires when its control bit is 1, which is exactly the behavior we need.
3. Repeat for the second qubit, the third, and so on.
4. Finally, reverse the order of the qubits with $\left\lfloor \frac{q}{2} \right\rfloor$ swap gates, since the product representation comes out with the bits reversed.

Counting gates: $q$ Hadamards, at most $q$ controlled rotations per qubit, and a handful of swaps, for a total of $O \left( q^2 \right)$ gates. Compare that to the FFT's $O \left( N \log N \right) = O \left( 2^q \cdot q \right)$: we went from exponential to quadratic in the number of bits. That is an exponential speedup for the transform itself.  
As a side note, the rotations $R_m$ for large $m$ are tiny (a phase of $\frac{2 \pi}{2^m}$), and in practice one can simply drop them with negligible effect on the result. This is called the *approximate QFT*, and it brings the gate count down to roughly $O \left( q \log q \right)$.

### The important caveat
If this sounds too good to be true, it partially is, and it is important to understand exactly why, since this is the trap most explanations skip.  
We do **not** get the Fourier coefficients. Amplitudes are not readable: measuring the register collapses it and hands us a single value $k$, drawn with probability $\left| y_k \right|^2$. We cannot look at the spectrum, we can only sample from it. We can of course run the entire thing again and get another sample, but reconstructing $N$ amplitudes from samples would take exponentially many runs, killing the whole advantage.

Luckily, sampling is precisely what period finding needs. Recall from the previous section that a periodic input makes the spectrum concentrate around the multiples of $\frac{N}{r}$, so a single measurement lands on (or very near) one of those spikes with a constant probability. Given such a $k$, the value $\frac{k}{N}$ is a good approximation of some $\frac{c}{r}$, and from that approximation we will be able to recover $r$. Let us now see how all of this fits together.

## Shor's algorithm
We now have all the ingredients. Let us assume we have already handled the easy classical cases (even $n$, prime $n$, prime powers), picked a random $a$ with $gcd \left( a, n \right) = 1$, and we wish to find the order $r$ of $a \pmod{n}$.

### The registers
We use two Quantum registers:
- The **first register** has $q$ qubits, where $q$ is chosen such that $n^2 \le N = 2^q \lt 2n^2$. This is the register that will hold $x$ and that we will apply the QFT to. The reason for the $n^2$ will become clear when we discuss continued fractions.
- The **second register** has $\left\lceil \log_2 n \right\rceil$ qubits, just enough to hold numbers modulo $n$. It will hold the values $a^x \pmod{n}$.

We start with the first register at $\left| 0 \right\rangle$ and the second register at $\left| 1 \right\rangle$ (the reason for starting at 1 rather than 0 will be clear in a moment).

### Step 1: superposition
We apply a Hadamard gate to each of the $q$ qubits of the first register, exactly like we did in Grover's algorithm. This creates the uniform superposition over all possible values of $x$:

$$\frac{1}{\sqrt{N}} \sum_{x=0}^{N-1} \left| x \right\rangle \left| 1 \right\rangle$$

### Step 2: modular exponentiation
Now we apply a unitary $U$ that does the following: $U \left| x \right\rangle \left| y \right\rangle = \left| x \right\rangle \left| y \cdot a^x \bmod n \right\rangle$. Since the second register starts at $y = 1$, the result is:

$$\frac{1}{\sqrt{N}} \sum_{x=0}^{N-1} \left| x \right\rangle \left| a^x \bmod n \right\rangle$$

Note $U$ is a legal Quantum operation: since $a$ is invertible modulo $n$, multiplying by $a^x$ is a permutation of the values modulo $n$, and permutations are reversible.  
How do we build it efficiently? With the same trick used in classical [modular exponentiation](https://en.wikipedia.org/wiki/Modular_exponentiation) (repeated squaring). Write $x$ in binary as $x = \sum_{l} x_l 2^l$, and then:

$$a^x = \prod_{l} \left( a^{2^l} \right)^{x_l}$$

The values $a^{2^l} \bmod n$ can be computed *classically* in advance (by squaring $q$ times). Then, for every qubit $x_l$ of the first register, we apply a "multiply the second register by the constant $a^{2^l} \bmod n$" circuit, controlled by $x_l$. That is $q$ controlled modular multiplications, each costing $O \left( \log^2 n \right)$ gates with the schoolbook method, for a total of $O \left( \log^3 n \right)$ gates.  
This is an important (and often hidden) point: the expensive part of Shor's algorithm is **not** the QFT. It is this modular exponentiation circuit, which is just classical arithmetic implemented reversibly.

It is tempting to say that we have now "computed $a^x \bmod n$ for all values of $x$ at once". That is true in a very limited sense: all the values are in there, but just like with the QFT, we cannot read them. If we measured now we would get a single random pair $\left( x, a^x \bmod n \right)$, which is something we could have easily computed classically. The value of this state is in its *structure*, not in its contents.

### Step 3: measuring the second register
Let us now measure the second register, and say we got some value $z$. Let $x_0$ be the smallest $x$ that yields $a^{x_0} \equiv z \pmod{n}$. Which values of $x$ are consistent with that measurement? Since the function has period $r$, exactly $x_0, x_0 + r, x_0 + 2r, \ldots$, and therefore the first register collapses to:

$$\frac{1}{\sqrt{M}} \sum_{j=0}^{M-1} \left| x_0 + jr \right\rangle$$

where $M \approx \frac{N}{r}$ is the number of such values below $N$.  
Look at that: the first register now holds a state that is **periodic with period $r$**, which is exactly the kind of input we said the Fourier transform is good at. The annoying part is the random offset $x_0$. We have no control over it, and if we measured the first register directly we would get $x_0 + jr$ for a random $j$, which tells us nothing about $r$. This is where the QFT comes in.  
(A technical remark: we do not actually have to measure the second register. Simply never touching it again gives the same statistics for the first register, a fact known as the [principle of deferred measurement](https://en.wikipedia.org/wiki/Deferred_Measurement_Principle). Measuring it just makes the explanation easier.)

### Step 4: QFT
We apply the QFT to the first register. Using the definition, the amplitude of $\left| k \right\rangle$ becomes:

$$y_k = \frac{1}{\sqrt{MN}} \sum_{j=0}^{M-1} e^{\frac{2 \pi i \left( x_0 + jr \right) k}{N}} = \frac{e^{\frac{2 \pi i x_0 k}{N}}}{\sqrt{MN}} \sum_{j=0}^{M-1} e^{\frac{2 \pi i j r k}{N}}$$

Here is the beautiful part: the offset $x_0$ moved out of the sum, and became a phase with absolute value 1. Since the probability of measuring $k$ is $\left| y_k \right|^2$, the offset **completely disappears** from the probabilities:

$$P \left( k \right) = \frac{1}{MN} \left| \sum_{j=0}^{M-1} e^{\frac{2 \pi i j r k}{N}} \right|^2$$

The QFT turned a shift in the value domain into a phase in the frequency domain, and phases are invisible to measurement. That is exactly why we needed it.

Let us look at the easy case first, when $r$ divides $N$. Then $M = \frac{N}{r}$ exactly, and the sum is the same roots-of-unity sum from before: it equals $M$ when $k$ is a multiple of $\frac{N}{r}$, and $0$ otherwise. Therefore, $P \left( k \right) = \frac{M^2}{MN} = \frac{1}{r}$ for each of the $r$ values $k = c \cdot \frac{N}{r}$ (with $c = 0, 1, \ldots, r-1$), and zero for everything else.

In the general case, $\frac{cN}{r}$ is not an integer, and the probability gets smeared around it. One can show (with some careful but unexciting bounding of the geometric series) that for every $c$, the integer $k$ closest to $\frac{cN}{r}$ is measured with probability of at least $\frac{4}{\pi^2 r}$. Summing over the $r$ possible values of $c$, we get that with probability of at least $\frac{4}{\pi^2} \approx 0.405$ we measure some $k$ with:

$$\left| k - \frac{cN}{r} \right| \le \frac{1}{2} \quad \Longrightarrow \quad \left| \frac{k}{N} - \frac{c}{r} \right| \le \frac{1}{2N}$$

### Step 5: measure
We measure the first register and get $k$. The Quantum part is now done, and everything from here on is classical. We know $k$ and $N$, and we know that $\frac{k}{N}$ is very close to some fraction $\frac{c}{r}$ with an unknown $c$ and an unknown $r \lt n$. How do we extract $r$?

## Continued fractions
A [continued fraction](https://en.wikipedia.org/wiki/Continued_fraction) is a way of writing a number as:

$$x = a_0 + \frac{1}{a_1 + \frac{1}{a_2 + \frac{1}{a_3 + \ldots}}}$$

which we write shortly as $\left[ a_0; a_1, a_2, a_3, \ldots \right]$. For a rational number like $\frac{k}{N}$ the expansion is finite, and computing it is nothing more than the [Euclidean algorithm](https://en.wikipedia.org/wiki/Euclidean_algorithm): take the integer part, invert the remainder, repeat. Cutting the expansion after each step gives a sequence of fractions called **convergents**, which are increasingly good approximations of the original number.

The reason continued fractions are useful for us is a classical theorem by Legendre: if $\left| x - \frac{c}{r} \right| \lt \frac{1}{2r^2}$, then $\frac{c}{r}$ (in lowest terms) is one of the convergents of $x$.  
Let us check we satisfy the condition. We know $\left| \frac{k}{N} - \frac{c}{r} \right| \le \frac{1}{2N}$, and we chose $N \ge n^2 \gt r^2$, so $\frac{1}{2N} \lt \frac{1}{2r^2}$. That is exactly why the first register needs about $2 \log_2 n$ qubits: it makes the approximation fine enough to be unambiguous (two different fractions with denominators smaller than $n$ are at least $\frac{1}{n^2}$ apart, so there is only one candidate close enough).

So, the classical recipe is: expand $\frac{k}{N}$ as a continued fraction, compute its convergents, and take the one with the largest denominator that is still smaller than $n$. That denominator is our candidate for $r$, which we verify by checking $a^r \equiv 1 \pmod{n}$ classically.

Let us see an example. Take $n = 21$ and $a = 2$. The order of 2 modulo 21 is 6 (the powers are $2, 4, 8, 16, 11, 1$), but we pretend we do not know that. Since $441 \le 512 \lt 882$ we use $q = 9$ qubits, so $N = 512$. Note 6 does not divide 512, so the peaks are around $\frac{512c}{6} \approx 0, 85.3, 170.7, 256, 341.3, 426.7$.  
Say we measured $k = 427$. Running the Euclidean algorithm:
- $427 = 0 \cdot 512 + 427$, so $a_0 = 0$.
- $512 = 1 \cdot 427 + 85$, so $a_1 = 1$.
- $427 = 5 \cdot 85 + 2$, so $a_2 = 5$.
- $85 = 42 \cdot 2 + 1$, so $a_3 = 42$.
- $2 = 2 \cdot 1 + 0$, so $a_4 = 2$.

Therefore $\frac{427}{512} = \left[ 0; 1, 5, 42, 2 \right]$, and the convergents are $\frac{0}{1}, \frac{1}{1}, \frac{5}{6}, \frac{211}{253}, \frac{427}{512}$. The largest denominator below 21 is 6, and indeed $2^6 = 64 \equiv 1 \pmod{21}$. We found $r = 6$.  
To finish the job: $r$ is even, and $2^3 = 8 \not\equiv -1 \pmod{21}$, so $gcd \left( 8 - 1, 21 \right) = 7$ and $gcd \left( 8 + 1, 21 \right) = 3$, and indeed $21 = 3 \cdot 7$.

### When things go wrong
There are two ways in which this step can fail:
1. **We measured $k$ far from any peak**, or measured $k = 0$ (which corresponds to $c = 0$ and tells us nothing). The verification $a^r \equiv 1 \pmod{n}$ catches the first case, and we simply run the Quantum part again.
2. **$c$ and $r$ share a common factor.** Continued fractions return $\frac{c}{r}$ *in lowest terms*, so if, say, $c = 2$ and $r = 6$, we get $\frac{1}{3}$ and a candidate $r' = 3$, which is a divisor of the real $r$. The verification fails since $a^{r'} \not\equiv 1 \pmod{n}$, but we got something useful.

This is where the promise from earlier comes in: we do not need to throw such results away. Every run gives us some divisor of $r$, so we can run the Quantum part a few times *with the same $a$* and take the [least common multiple](https://en.wikipedia.org/wiki/Least_common_multiple) of the denominators we got. Two runs already yield the correct $r$ with probability of at least $\frac{1}{4}$, and every additional run only improves the odds. Another cheap trick is trying small multiples of a candidate $r'$ (e.g. $2r'$, $3r'$) before running the Quantum part again.

## A complete example: factoring 15
Let us put it all together on the smallest interesting number, $n = 15$. It is odd, it is not prime and it is not a prime power, so we proceed.
1. We pick a random $a$, say $a = 7$. Since $gcd \left( 7, 15 \right) = 1$, we did not get lucky, and we need to find the order.
2. Since $225 \le 256 \lt 450$ we use $q = 8$ qubits for the first register, so $N = 256$. The second register needs 4 qubits to hold numbers modulo 15.
3. After the Hadamards and the modular exponentiation we have $\frac{1}{16} \sum_{x=0}^{255} \left| x \right\rangle \left| 7^x \bmod 15 \right\rangle$. Since $7^0, 7^1, 7^2, 7^3, 7^4 \equiv 1, 7, 4, 13, 1 \pmod{15}$, we can group the terms by the value of the second register:

$$\frac{1}{16} \left[ \left( \left| 0 \right\rangle + \left| 4 \right\rangle + \ldots + \left| 252 \right\rangle \right) \left| 1 \right\rangle + \left( \left| 1 \right\rangle + \left| 5 \right\rangle + \ldots + \left| 253 \right\rangle \right) \left| 7 \right\rangle + \left( \left| 2 \right\rangle + \left| 6 \right\rangle + \ldots + \left| 254 \right\rangle \right) \left| 4 \right\rangle + \left( \left| 3 \right\rangle + \left| 7 \right\rangle + \ldots + \left| 255 \right\rangle \right) \left| 13 \right\rangle \right]$$

4. We measure the second register, and say we got 4 (each of the 4 values is equally likely). The first register collapses to $\frac{1}{8} \sum_{j=0}^{63} \left| 2 + 4j \right\rangle$, which is periodic with period 4 and offset 2.
5. We apply the QFT. Since $r = 4$ divides $N = 256$ we are in the easy case: the amplitude of $\left| k \right\rangle$ is $\frac{1}{128} e^{\frac{2 \pi i \cdot 2k}{256}} \sum_{j=0}^{63} e^{\frac{2 \pi i j k}{64}}$, which is zero unless $k$ is a multiple of 64, in which case it has absolute value $\frac{64}{128} = \frac{1}{2}$. So we measure each of $k = 0, 64, 128, 192$ with probability $\frac{1}{4}$, and notice the offset 2 had no effect on that whatsoever.
6. Let us go over the possible outcomes:
   - $k = 0$: gives $\frac{0}{1}$, which is useless. We run again.
   - $k = 64$: gives $\frac{64}{256} = \frac{1}{4}$, so $r = 4$. Verified, since $7^4 \equiv 1 \pmod{15}$.
   - $k = 128$: gives $\frac{1}{2}$, so $r' = 2$. This is the case of $c = 2$ sharing a factor with $r = 4$. Verification fails ($7^2 \equiv 4 \pmod{15}$), but trying $2r' = 4$ succeeds.
   - $k = 192$: gives $\frac{3}{4}$, so $r = 4$. Verified.
7. With $r = 4$: it is even, and $7^2 = 49 \equiv 4 \pmod{15}$, which is not $-1 \equiv 14 \pmod{15}$. Therefore we get $gcd \left( 4 - 1, 15 \right) = 3$ and $gcd \left( 4 + 1, 15 \right) = 5$, and indeed $15 = 3 \cdot 5$.

## Summary and complexity
Let us write down Shor's algorithm in full:
1. If $n$ is even, prime or a prime power, handle it classically.
2. Pick a random $1 \lt a \lt n - 1$. If $gcd \left( a, n \right) \gt 1$, we found a factor.
3. Use the Quantum computer to find the order $r$ of $a \pmod{n}$:
   - Prepare $\frac{1}{\sqrt{N}} \sum_x \left| x \right\rangle \left| 1 \right\rangle$ with $n^2 \le N \lt 2n^2$.
   - Apply modular exponentiation to get $\frac{1}{\sqrt{N}} \sum_x \left| x \right\rangle \left| a^x \bmod n \right\rangle$.
   - Apply the QFT to the first register and measure it to get $k$.
   - Use continued fractions on $\frac{k}{N}$ to get a candidate $r$, verify it, and repeat (combining results with LCM) if needed.
4. If $r$ is odd or $a^{\frac{r}{2}} \equiv -1 \pmod{n}$, go back to step 2.
5. Otherwise, $gcd \left( a^{\frac{r}{2}} \pm 1, n \right)$ are non-trivial factors of $n$.

To summarize the whole trick in words: we use superposition and modular exponentiation to create a state that is periodic in $x$ with period $r$, but with a random offset we cannot control; we apply the QFT to convert that period into sharp frequencies, turning the offset into an invisible phase; and we measure once to sample one of those frequencies. The Quantum computer never tells us the entire Fourier spectrum. It just tells us the one thing we actually wanted to know.

How fast is it? All the steps are polynomial in the number of bits of $n$: the modular exponentiation dominates with $O \left( \log^3 n \right)$ gates (less with fancier multiplication), the QFT needs only $O \left( \log^2 n \right)$, and the expected number of repetitions is a constant. Compare that to the best known classical algorithm, the [General Number Field Sieve](https://en.wikipedia.org/wiki/General_number_field_sieve), whose running time is roughly $e^{1.9 \left( \ln n \right)^{\frac{1}{3}} \left( \ln \ln n \right)^{\frac{2}{3}}}$, which is super-polynomial in the number of bits.

## What it means for RSA
Since RSA's security relies on factoring being hard, Shor's algorithm breaks RSA, at least in theory. In practice, at the time of writing, the largest numbers "factored" by running Shor's algorithm on real Quantum hardware are tiny (15 and 21), and even those demonstrations used simplified circuits that relied on knowing the answer in advance.  
The problem is noise: the circuits we described need many qubits that stay coherent for billions of operations, which requires [Quantum error correction](https://en.wikipedia.org/wiki/Quantum_error_correction) and therefore many physical qubits per logical qubit. Resource estimates have been dropping steadily, though. A 2025 estimate by Craig Gidney puts factoring a 2048-bit RSA modulus at under a million noisy qubits running for less than a week. That is still far beyond today's machines, but no longer absurd.  
This is why the cryptographic community has not been waiting: NIST has already standardized [post-quantum algorithms](https://en.wikipedia.org/wiki/Post-quantum_cryptography) (such as ML-KEM), which are based on problems not known to be vulnerable to Shor's algorithm, and migration to them is already underway.

Stay tuned!

Jonathan Bar Or (https://jonathanbaror.com)
