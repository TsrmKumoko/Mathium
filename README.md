# Mathium

Mathium is a python library which supports complex numbers, fractions, matrices, polynomials calculation.

## Intro

There are 6 main classes in Mathium: Real, Complex, Matrix, Tensor, Fraction, Polynomial. They can be freely mixed up when you calculate a certain expression.

Use the method `display` to view the output, which will be rendered by $\LaTeX$. Use the method `LaTeX` to get the LaTeX code of the output expression.

Basic methods:

| Operator | Real | Complex | Matrix | Tensor | Fraction | Polynomial |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `a + b` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `a - b` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `a * b` | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ |
| ` a / b` | ✅ | ✅ | ❌ | ❌ | ❌ | ✅ |
| `a ** b` | ❌ | ❌ | ✅ | ❌ | ❌ | ❌ |
| `a[n]` | ❌ | ❌ | ❌ | ✅ | ❌ | ❌ |

Special methods:

| Real | Complex | Matrix | Tensor | Fraction | Polynomial |
| :---: | :---: | :---: | :---: | :---: | :---: |
| `factors(self)` | `Re(self)`<br>`Im(self)`<br>`conj(self)`<br>`norm(self)` | `I(int)`<br>`det(Matrix)`<br>`adjoint(Matrix)`<br>`inv(Matrix)` | `LineToCoor(self, idx)`<br>`PolyConv(ten0, ten1, map0, map1)`<br>`PolyPlus(ten0, ten1, map0, map1)` | `reduce(self)` | |

## Examples

### Real

The `Real` class is used to unify the calculation, so it's not so meaningful to directly use it. 

You can factorize a integer by using `Real.factors()`:

```python
>>> Real(36).factors()
[2, 2, 3, 3]
```

### Complex

Give the real and imaginary part of a complex number to initialize a `Complex` object.

```python
c = Complex(1, 1)
display(c * Complex(1, -2))
```

You can get a LaTeX object which should be viewed with Jupyter. That is:

$$
3-\mathrm{i}
$$

If you want to get the LaTeX code, do:
```python
>>> LaTeX(Complex(1, 1))
'1-\\mathrm{i}'
```

### Matrix

```python
m1 = Matrix(
    [1, 2],
    [-1, Complex(1, 1)]
)
m2 = Matrix(
    [1, Complex(1, 1)],
    [3, Complex(1, -1)]
)
display(m1 * m2)
```
$$
\begin{bmatrix}
7 & 3-\mathrm{i} \\
2+3\mathrm{i} & 1-\mathrm{i}
\end{bmatrix}
$$

### Tensor

This class is *temporarily* made for `Polynomial` class. Using it directly is not suggested now because many features have not yet been implemented.

```python
t = Tensor([
    [0, 1],
    [2, 3]
])
t[0][1]
```
$$
1
$$

### Fraction

Now float can be casted to fraction. Reduction is available and will be automatically done. Use `auto_reduce = False` to turn it off.

```python
r1 = Fraction(0.128)
r2 = Fraction(25, 7)
display(r1 * r2)
```
$$
\frac{16}{35}
$$

### Polynomial

```python
t = Polynomial([1, 1], 'w') / Polynomial([-1, 1], 'w')
display(t*t - t + 2)

f1 = Polynomial([
    [1, 1],
    [1, 1]
], ['x', 'y'])
f2 = Polynomial(
    [1, -1, 1],
    ['x']
)
display(f1 * f2)
```
$$
\frac{2w^{3}-4w^{2}+6w-4}{w^{3}-3w^{2}+3w-1}
$$

$$
x^{3}y+x^{3}+y+1
$$

## Future Version

There are a number of features undone which will be finished little by little. Updates will include:

* Factorization in `Polynomial` class
* Improved `Matrix` and `Tensor` class
* Derivative
