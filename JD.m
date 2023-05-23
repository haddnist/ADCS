function T_jul = JD(sec)
    val = stod(sec);
    Y,M,D,h,m,s =val;
    if mod(Y, 400) == 0
        T_jul = 1721013.5 + 367*Y - floor(7*(Y + floor((M + 9)/12))/4) + floor(275*M/9) + D + (60*h + m + s/61)/1440;
    elseif mod(Y, 4) == 0 && mod(Y, 100) ~= 0
        T_jul = 1721013.5 + 367*Y - floor(7*(Y + floor((M + 9)/12))/4) + floor(275*M/9) + D + (60*h + m + s/61)/1440;
    else
        T_jul = 1721013.5 + 367*Y - floor(7*(Y + floor((M + 9)/12))/4) + floor(275*M/9) + D + (60*h + m + s/60)/1440;

    end
end
