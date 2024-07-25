module Rotations

using LinearAlgebra

function rotate(angle::Float64, axis::Union{Char, Int})
    """
    Coordinate rotation about a principal axis.
    :param angle: magnitude of angle of rotation (radians)
    :param axis: axis to rotate about: 'x', 'y', 'z'; or 1, 2, 3
    :return matrix: rotation matrix of "angle" about an axis 
    """
    if axis == 'x' || axis == 1
        matrix = [
            1.0  0.0           0.0;
            0.0  cos(angle)  sin(angle);
            0.0 -sin(angle)  cos(angle)
        ]
    elseif axis == 'y' || axis == 2
        matrix = [
            cos(angle)  0.0 -sin(angle);
            0.0         1.0  0.0;
            sin(angle)  0.0  cos(angle)
        ]
    elseif axis == 'z' || axis == 3
        matrix = [
            cos(angle)  sin(angle)  0.0;
           -sin(angle)  cos(angle)  0.0;
            0.0         0.0         1.0
        ]
    else
        throw(ArgumentError("Invalid axis"))
    end
    return matrix
end

end # module
