include ("derivative.jl")

function quadratic(f)
        #compute first derivative of f
        f1=derivative (f)

        #get y intercept
        c=f(0.0)

        #get y intercept of first derivative
        b=f1(0.0)

        a=f(1.0)-b-c

        return (-b+sqrt(b^2 -4a*c + 0im))/2a, (-b-sqrt(b^2 -4a*c +0im))/2a
        end

