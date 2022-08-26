using ojwul
using StaticArrays
using Test

function testcamera(cam::AbstractCamera)
    x = @SVector randn(2)
    err = (@SVector randn(2)) * 1.e-6
    xe = x .+ err
    W = @SMatrix randn(2, 2)

    # Test image to ideal transformations
    y, Wy = image2ideal(cam, x, W)
    @test image2ideal(cam, x) == y
    @test isapprox(ideal2image(cam, y), x)

    # Test warping of the weight matrix
    ye = image2ideal(cam, xe)
    @test isapprox(Wy * (ye - y), W * err; rtol=1.e-6)
end

@testset "camera.jl" begin
    halfimsz = SA[640, 480]
    x = @SVector randn(2)
    x = x .* (0.3 * halfimsz)
    err = (@SVector randn(2)) * 1.e-6
    xe = x .+ err
    W = @SMatrix randn(2, 2)

    # Test pixel to image transformations
    y, Wy = pixel2image(halfimsz, x, W)
    @test pixel2image(halfimsz, x) == y
    @test isapprox(image2pixel(halfimsz, y), x)

    # Test warping of the weight matrix
    ye = pixel2image(halfimsz, xe)
    @test isapprox(Wy * (ye - y), W * err; rtol=1.e-6)

    # Test cameras
    testcamera(SimpleCamera(abs(randn())))
    f = abs.(@SVector randn(2))
    c = @SVector randn(2)
    testcamera(NoDistortionCamera(f, c))
    testcamera(ExtendedUnifiedCamera(f, c, rand(), abs(randn()) * 0.1))
end
