# Mathium

Mathium is a python library which supports complex numbers, fractions, matrices, polynomials calculation.

## Intro

There are 5 main classes in Mathium: Complex, Matrix, Tensor, Fraction, Polynomial. They can be freely mixed up when you calculate a certain expression.

Use the method `display` to view the output, which will be rendered by $\LaTeX$. Use the method `LaTeX` to get the LaTeX code of the output expression.

Basic methods:

| Operator | Complex | Matrix | Tensor | Fraction | Polynomial |
| :---: | :---: | :---: | :---: | :---: | :---: |
| `a + b` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `a - b` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `a * b` | ✅ | ✅ | ❌ | ✅ | ✅ |
| ` a / b` | ✅ | ❌ | ❌ | ❌ | ✅ |
| `a ** b` | ❌ | ✅ | ❌ | ❌ | ❌ |
| `a[n]` | ❌ | ❌ | ✅ | ❌ | ❌ |

Special methods:

| Complex | Matrix | Tensor | Fraction | Polynomial |
| :---: | :---: | :---: | :---: | :---: |
| `Re(self)`<br>`Im(self)`<br>`conj(self)`<br>`norm(self)` | `I(int)`<br>`det(Matrix)`<br>`adjoint(Matrix)`<br>`inv(Matrix)` | `LineToCoor(self, idx)`<br>`PolyConv(ten0, ten1, map0, map1)`<br>`PolyPlus(ten0, ten1, map0, map1)` | `reduce(self)` | |

## Usage

### Complex

Give the real and imaginary part of a complex number to initialize a `Complex` object.

```python
c = Complex(1, 1)
display(c * Complex(1, -2))
```

You can get a LaTeX object. Please view it by Jupyter. That is:
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
Output:
```python
1
```

### Fraction

```python
r1 = Fraction(1, 2)
r2 = Fraction(2, 3)
display(r1 + r2)
```
$$
\frac{7}{6}
$$

### Polynomial

```python
t = Polynomial([1, 1], 'w') / Polynomial([-1, 1], 'w')
display(t*t - t + 2)

f1 = Polynomial([
    [3, 2],
    [1, -1]
], ['x', 'y'])
f2 = Polynomial([3, 0, 1], ['z'])
display(f1 * (2 * f2 - f1))
```
$$
\frac{2w^{3}-4w^{2}+6w-4}{w^{3}-3w^{2}+3w-1}
$$
$$
-x^{2}y^{2}+2x^{2}y-x^{2}+4xy^{2}-2xyz^{2}-4xy+2xz^{2}-4y^{2}+4yz^{2}+6z^{2}+9
$$

There are a number of features undone which will be finished little by little.
