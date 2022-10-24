function print_to_serial(serial_dev, st, use_fprintf)

disp(st) ; 

if use_fprintf
    %\r = enter
    st = [st, '\r'] ;
    fprintf(serial_dev,st);
else
    for i=1:length(st)
        fwrite(serial_dev,st(i),'uchar') ; 
    end
end

end
