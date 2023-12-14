function bool = varA_check(A)

    if isempty(find(abs(A)<1e-18,1))
        bool = 1;
    else 
        bool = 0;
    end

end



    