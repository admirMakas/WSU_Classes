L = length(Time_chan_2);
j=0;
for i=1:L
    if i>1
        if Time_chan_2(i-1)<Time_chan_2(i) &&  Time_chan_2(i+1)<Time_chan_2(i)
            j=j+1;
            MaxVal1(j, 1) = Time_chan_2(i);
            MaxVal1(j, 2) = Time_domain(i);
        end
    end
end
j=0;
N = length(MaxVal1);
for i=1:N
    if i>1
        if MaxVal1(i-1)<MaxVal1(i) &&  MaxVal1(i+1)<MaxVal1(i)
            j=j+1;
            MaxVal2(j, 1) = MaxVal1(i, 1);
            MaxVal2(j, 2) = MaxVal1(i, 2);
        end
    end
end