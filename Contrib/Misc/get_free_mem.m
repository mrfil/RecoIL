function [mem, unit] =get_free_mem()
%Usage: [mem, unit] =get_free_mem()

    [~,out]=system('vmstat -s -S M | grep "free memory"');
    
    mem=sscanf(out,'%f  free memory');

    unit = 'GB' ; 
    
    mem = mem/1000; % Convert to GB.

end