
function QACEval(tn, q, qd, SMDT, STSDAT, par, p_contact)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    flag_contact,damp=p_contact
    uz = [0; 0; 1]
    uy = [0; 1; 0]

    QAC = zeros(ngc)
    if flag_contact == 1
        # Account for gravitational force in negative y direction
        for i in [1,]# body number
            # mi = SMDT[1, i]
            QACi = vcat(-qd[1]* damp, -qd[2] * damp, -qd[3]* damp, zeros(4))
            QAC = add_constraint!(QAC, QACi, 7 * (i - 1), 0)
        end
    end
    if q[3]<-0.5
        println("q[3] is less than -0.5,",(q[3]+0.5))
        QACi = vcat(0, 0, 1e9*(abs(q[3]+0.5))^1.5, zeros(4))
        QAC = add_constraint!(QAC, QACi, 7 * (1 - 1), 0)
    end
    return QAC
end