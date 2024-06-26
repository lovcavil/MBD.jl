export RevPart,TranPart,RotDrPart,CylPart,DistDrPart,DistPart,SphPart
function RevPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    vx1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    vz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    vx2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]
    vz2pr = [SJDT[20, k], SJDT[21, k], SJDT[22, k]]

    return i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr
end
function TranPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    ux1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    uz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    ux2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]
    uz2pr = [SJDT[20, k], SJDT[21, k], SJDT[22, k]]

    return i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr
end
function RotDrPart(k, SJDT)
    # Data derived from the host revolute or cylindrical joint data table
    i = SJDT[2, k]
    j = SJDT[3, k]
    vx1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    vz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    vy1pr = cross(vz1pr, vx1pr)  # using cross product to get the perpendicular vector
    vx2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]

    return i, j, vx1pr, vy1pr, vx2pr
end
function CylPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    ux1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    uz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    ux2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]
    uz2pr = [SJDT[20, k], SJDT[21, k], SJDT[22, k]]

    return i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr
end
function DistDrPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    return i, j, s1pr, s2pr
end
function DistPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    d = SJDT[10, k]
    return i, j, s1pr, s2pr, d
end
function SphPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]

    return i, j, s1pr, s2pr
end