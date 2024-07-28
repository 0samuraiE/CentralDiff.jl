using CentralDiff
using Test

@testset "CentralDiff.jl" begin
    _centered_coords(x0, h, N) = @. x0 + h * ((1-N)//2:(N-1)//2)

    for O in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
        f̃(x) = exp(x)
        df̃dx(x) = exp(x)
        d2f̃dx(x) = exp(x)
        x0 = 1.0

        order = CentralDiff.countorder(2, 1, 1 / 2) do h
            X1 = _centered_coords(x0, h, O)
            Y1 = f̃.(X1)

            X2 = _centered_coords(x0, h, 2O - 1)
            Y2 = f̃.(X2)

            (;
                fc=f̃(x0) - fc(Order(O), Y1),
                dfdxc=df̃dx(x0) - dfdxc(Order(O), Y1, 1 / h),
                dfdx=df̃dx(x0) - dfdx(Order(O), Y2, 1 / h),
                d2fdx=d2f̃dx(x0) - d2fdx(Order(O), Y2, 1 / h),
            )
        end
        @test all(order.fc .> O)
        @test all(order.dfdxc .> O)
        @test all(order.dfdx .> O)
        @test all(order.d2fdx .> O)
    end

    for O in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
        f̃(x) = exp(x)
        df̃dx(x) = exp(x)
        d2f̃dx(x) = exp(x)
        x0 = 1.0

        order = CentralDiff.countorder(2, 1, 1 / 2) do h
            X = _centered_coords(x0, h, O + 1)
            Y = f̃.(X)

            (;#
                dfdx=df̃dx(x0) - Simple.dfdx(Order(O), Y, 1 / h),
                d2fdx=d2f̃dx(x0) - Simple.d2fdx(Order(O), Y, 1 / h),
            )
        end
        @test all(order.dfdx .> O)
        @test all(order.d2fdx .> O)
    end
end
