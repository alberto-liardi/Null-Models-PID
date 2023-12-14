function [] = PIDcheck(PID)
    % values above -1e-11 are considered non-negative

    % in case of a matric of PID atoms (each row consists of PID atoms for
    % a specific redundancy function)
    for s = 1:size(PID,1)

        Ux = PID(s,1);
        Uy = PID(s,2);
        R = PID(s,3);
        S = PID(s,4);
        
        MI = Ux + Uy + R + S;
        MIx = Ux + R;
        MIy = Uy + R;
        assert(MIx > -1e-11 && MIy > -1e-11, "Single mutual information is negative");
        assert(MI > max(MIx, MIy) - 5e-10, "MI is smaller than a singular MI");
        assert(S > -1e-11, "Synergy is negative!");

    end
end