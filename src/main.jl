using Intervals
using Plots

# gauss-legendre (2 points)

function int(f, intersection)
    a = first(intersection)
    b = last(intersection)
    x = 1 / sqrt(3)
    if b - a != 0
        (
            f(((b - a) / 2) * (x) + ((a + b) / 2)) +
            f(((b - a) / 2) * (-x) + ((a + b) / 2))
        ) * ((b - a) / 2)
    else
        0
    end
end


# base functions

function e_base(i, x, h)
    x_a = h * (i - 1)
    x_b = h * i
    x_c = h * (i + 1)

    if x_a <= x <= x_b
        (x - x_a) / h
    elseif x_b <= x <= x_c
        (x_c - x) / h
    else
        0
    end
end


function e_deriv(i, x, h)
    x_a = h * (i - 1)
    x_b = h * i
    x_c = h * (i + 1)

    if x_a <= x <= x_b
        1 / h
    elseif x_b <= x <= x_c
        -1 / h
    else
        0
    end
end

# B (left side)
B(w, dw, v, dv, a, b) = w(0) * v(0) - int(x -> dw(x) * dv(x), intersect(a .. b, 0 .. 3))

# L (right side)
L(v, a, b) =
    -1 / 10 * int(v, intersect(a .. b, 0 .. 1)) -
    1 / 5 * int(v, intersect(a .. b, 1 .. 2)) - int(v, intersect(a .. b, 2 .. 3)) + 5v(0)


function main(n)
    
    a = 0
    b = 3

    h = (b - a) / n

    e(i) = x -> e_base(i, x, h)
    de(i) = x -> e_deriv(i, x, h)

    shift = x -> 2 * e_base(n, x, h)
    dshift = x -> 2 * e_deriv(n, x, h)


    A_mat = zeros(n, n)

    # smart fill

    mainDiagonalValue = B(e(1), de(1), e(1), de(1), h * (1 - 1), h * (1 + 1))
    sideDiagonalValue = B(e(1), de(1), e(0), de(0), h * (1 - 1), h * (0 + 1))

    for i in range(1, n)
        A_mat[i, i] = mainDiagonalValue
    end

    for i in range(2, n)
        A_mat[i, i-1] = sideDiagonalValue
        A_mat[i-1, i] = sideDiagonalValue
    end

    A_mat[1, 1] = B(e(0), de(0), e(0), de(0), h * (0 - 1), h * (0 + 1))

    # calculate every

    # for i in range(1, n)
    #     A_mat[i, i] =
    #         B(e(i - 1), de(i - 1), e(i - 1), de(i - 1), h * (i - 1 - 1), h * (i - 1 + 1))
    # end

    # for i in range(2, n)
    #     A_mat[i, i-1] = B(
    #         e(i - 1),
    #         de(i - 1),
    #         e(i - 1 - 1),
    #         de(i - 1 - 1),
    #         h * (i - 1 - 1),
    #         h * ((i - 1 - 1) + 1)
    #     )
    #     A_mat[i-1, i]  = B(
    #         e(i - 1 - 1 ),
    #         de(i - 1- 1),
    #         e(i - 1 ),
    #         de(i - 1 ),
    #         h * (i - 1 - 1),
    #         h * ((i - 1 - 1) + 1))
    # end

    # A_mat[1, 1] = B(e(0), de(0), e(0), de(0), h * (0 - 1), h * (0 + 1))


    B_mat = zeros(n, 1)

    for i in range(1, n)
        B_mat[i, 1] = L(e(i - 1), h * (i - 1 - 1), h * (i - 1 + 1))
    end

    B_mat[n, 1] -= B(shift, dshift, e(n - 1), de(n - 1), h * (n - 1), h * (n - 1 + 1))

    U_mat = A_mat \ B_mat


    
    display(U_mat)
    display(A_mat)
    display(B_mat)

    u = x -> shift(x) + sum(U_mat[i, 1] * e(i - 1)(x) for i = 1:n)



    der2(f, h) = x -> (f(x+h)+f(x-h)-2f(x))/h^2
   
    println(der2(u, h)(0.5))
    println(der2(u, h)(1.5))
    println(der2(u, h)(2.5))

    display(plot(0.01:h:2.99, der2(u, h)))
    display(plot(0:h:3, u))
    nothing
end

main(500)