include("lgquad.jl")

macro run()
    :(include("main.jl") ; main())
end

function main()
    # for now there are performance issues with n > 30
    n = 4;
    f(x) = x^2
    a = 0.0
    b = 1.0

    println("integral: $(lgquad(f, n, a, b))")
end