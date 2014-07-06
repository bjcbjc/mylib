function survival = calSurvival(days2death, days2lastfollow, vitalstatus)

    n = length(days2death);
    survival = NaN(n, 2); %survival length, censor
    
    i = strcmpi(vitalstatus, 'DECEASED');
    survival(i, 1) = days2death(i);
    survival(i, 2) = 0;
    
    i = strcmpi(vitalstatus, 'LIVING');
    survival(i, 1) = days2lastfollow(i);
    survival(i, 2) = 1;

end



