function clear_cache()
%Clears the hard drive read cache used by load_image()
%
%   === Notes ===
%   This function is not strictly necessary but is included for
%   convenience.  The cache used by load_image() can be cleared by calling
%   that function and specifying the option 'clear_cache' to be true.  This
%   function simply calls load_image() with 'clear_cache' set to true by
%   providing it a fake filename to load, then swallowing the generated
%   error.
%
%   === Example Usage ===
%   >> clear_cache();


try
    %Clears cache. Also throws an error because '' is not a real filename
    load_image('',true);
catch
end

end