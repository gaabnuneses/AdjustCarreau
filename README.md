# AdjustCarreau
Fit a Carreau curve to a set of measured points in a flowcurve using Julia language.

## Initialization
Assign a column vector to each variable x,y.

```julia
using Plots # Visualize your data
x,y = Array{Float64,1},Array{Float64,1}

# x → usually the shear rate in [1/s]
# y → usually the dynamic viscosity in [Pa.s]
```

## Create some functions
1. Function to limit your data on the x axis

```julia
function Limite(x,y,xmin,xmax)
    X = []
    Y = []
    for i in 1:length(x)
        if x[i]>xmin && x[i]< xmax
            X=[X;x[i]]
            Y=[Y;y[i]]
        end
    end
    return X,Y #
end
```
2. Function to calculate error energy of descendant gradient method

```julia
function erro(x,y,n,λ,μ₀,μ∞)
    N = length(x)
    # Sum initilization for each error derivative. 
    ϵn = 0
    ϵλ = 0
    ϵμ₀ = 0
    ϵμ∞ = 0
    # (x,y) dependant terms of the sum. 
    for i in 1:N
        ϵn += ((μ∞-y[i])*(λ^2*x[i]^2+1)^((1/2)*n-1/2)+(μ₀-μ∞)*(λ^2*x[i]^2+1)^(n-1))*(μ₀-μ∞)*log(λ^2*x[i]^2+1)
        ϵλ += x[i]^2*(μ₀-μ∞)*((μ∞-y[i])*(λ^2*x[i]^2+1)^((1/2)*n-3/2)+(λ^2*x[i]^2+1)^(n-2)*(μ₀-μ∞))
        ϵμ₀ += (μ∞-y[i])*(λ^2*x[i]^2+1)^((1/2)*n-1/2)+(μ₀-μ∞)*(λ^2*x[i]^2+1)^(n-1)
        ϵμ∞ +=  (μ∞-y[i])*(λ^2*x[i]^2+1)^((1/2)*n-1/2)+(μ₀-μ∞)*(λ^2*x[i]^2+1)^(n-1)
    end
    
    # Independent terms addition
    ϵλ = ϵλ*(2*(n-1))*λ
    ϵμ₀ = ϵμ₀*2
    ϵμ∞ = ϵμ∞*2 + (2*N+2)*μ∞
    
    return ϵn,ϵλ,ϵμ₀,ϵμ∞
end
```

3. Coefficient correction function

```julia
function calCoef(x,y,n,λ,μ₀,μ∞,nit,h=.01)
    for i in 1:nit
        # Call the previous function
        ϵn,ϵλ,ϵμ₀,ϵμ∞ = erro(x,y,n,λ,μ₀,μ∞)
        # Correct the terms with underrelaxation coefficient 0>h>1
        n = n - h*ϵn
        λ = λ - h*ϵλ
        μ₀ = μ₀ - h*ϵμ₀
        μ∞ = μ∞ - h*ϵμ∞
    end
    return n,λ,μ₀,μ∞
end
```

4. Normalized adjust

```julia
function adjustCarreau(x,y)
    # Set data limits as you wish
    x,y=Limite(x,y,100,330)
    
    # Find normalization terms
    xmx =maximum(x)
    ymx = maximum(y)
    ymn = minimum(y)
    
    # x-axis normalization
    x_n = (x)./(maximum(x))
    # y-axis normalization
    y_n = (y .- minimum(y))./(maximum(y) -minimum(y))
    
    # Find normalized coefficients
    n,λ,μ₀,μ∞ = calCoef(x_n,y_n,-2,3,.75,0.25,30000,.0001)
    
    # Unormalize Carreau Coefficients
    M∞ = μ∞*(ymx-ymn)+ymn
    M₀ = μ₀*(ymx-ymn)+ymn
    L = λ/xmx
    N = n
    
    # Show on REPL
    println("viscosidade final μ∞ = $(M∞) Pa.s")
    println("viscosidade inicial μ₀ =$(M₀) Pa.s")
    println("Tempo de relaxação λ = $L s")
    println("Coeficiente power law n = $N")

    # Plot
    xx = 0:5:xmx
    yy(x) = M∞ + (M₀-M∞)*(1 +(L*x)^2)^((N-1)*(1/2))
    p = plot!(xx,yy.(xx),lw=6,label="Ajuste Carreau")
    return p
end
```

5. Do your adjust calling the function and save the results.

```julia
p = adjustCarreau(x,y)
savefig(p,"Adjust.png")

```
