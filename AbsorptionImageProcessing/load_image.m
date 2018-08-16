function image_out = load_image(filename, clear_cache)
%Loads data from the csv ascii files produced by our camera software
%   === Inputs ===
%   filename should be the name of the file with extension.  A relative or
%   absolute path may also be included in filename.
%   
%   clear_cache (optional) tells the function whether or not to clear data
%   from the hard drive read cache.  By default this is false.
%
%   === Outputs ===
%   image_out is a 2D array of data from the specified file.
%
%   === Notes ===
%   Results from reading the hard drive are cached for speed.  Call this
%   function with clear_cache=true to clear the cache.  This is necessary
%   if, for example, data on the hard drive has changed and the new file
%   should be read in.  Also, once the cache is full, no more files are
%   added to the cache (currently we don't delete old cache entries to add
%   new ones). If you'd like to add new files to the cache but it is full,
%   you must clear it.  Caching of the inputs is done by treating it as a
%   raw string, so if two different paths to the same file are specified in
%   two different calls, the function will cache the result twice, once for
%   each given input string.
%
%   === Example Usage ===
%   >> filename = fullfile('20170412','Savefile_45_back.ascii');
%   >> image_out = load_image(filename); %load an image
%   >> plot_image(image_out,'Test of load_image()');
%   >> %To clear the cache, run the following line (will error since a  
%   >> %nonexsistent file name is provided, but the cache will still be cleared)
%   >> image = load_image('', true); 

%Create persistent variable for storing memoized function
persistent cached_reader;
max_cache_size=10000; %Maximum number of files to cache

%If this is the first call, create the memoized function object
if isempty(cached_reader) 
    cached_reader=containers.Map;
end

%Clear the cache if instructed
if nargin<2
    clear_cache=false;
end
if clear_cache
    %Clear the memoization cache
    cached_reader=containers.Map;
    remove(cached_reader, keys(cached_reader) );
end

%Get the image data, using the cache if possible
if isKey(cached_reader,filename)
    %We have a cached result, so we'll return it
    image_out=cached_reader(filename);
else
    image_out=read_image_from_drive(filename);
    %add to cache if there's room
    if length(cached_reader)<max_cache_size
        cached_reader(filename)=image_out;
    end
end

%Uncomment the following line and comment the above code to ditch the
%memoizing stuff and just read from the hard drive
% image_out=read_image_from_drive(filename);

end

function image_out = read_image_from_drive(filename)
%Loads data from the files produced by our camera software
%   === Inputs ===
%   filename should be the name of the file with extension.  A relative or
%   absolute path may also be included in filename.
%
%   === Outputs ===
%   image_out is a 2D array of data from the specified file.
%
%   === Notes ===
%   This function is designed to load in tab-separated csv files that end
%   in '.ascii' or png files that end in '.png'.
%
%   === Example Usage ===
%   This is a helper function that cannot be called directly.

%Check if the file name ends in '.ascii', and if so call the appropriate
%subfunction
[~,~,extension]=fileparts(filename);
if strcmp(extension,'.ascii')
    image_out=read_image_from_drive_ascii(filename);
else
    image_out=read_image_from_drive_png(filename);
end

end

function image_out = read_image_from_drive_png(filename)
%Loads data from a png file produced by our camera software
%   === Inputs ===
%   filename should be the name of the file with extension.  A relative or
%   absolute path may also be included in filename.
%
%   === Outputs ===
%   image_out is a 2D array of data from the specified file.
%
%   === Notes ===
%   This function is designed only for grayscale png images, typically with
%   16 bit depth. It may work for other formats as well as long as imread
%   knows how to deal with them, but the images should still be gray scale
%
%   === Example Usage ===
%   This is a helper function that cannot be called directly.

%Call imread to read in the file
image_out=imread(filename);

%Convert data type to double precision, necessary for calling eigs
image_out=double(image_out);

end

function image_out = read_image_from_drive_ascii(filename)
%Loads data from the csv ascii files produced by our camera software
%   === Inputs ===
%   filename should be the name of the file with extension.  A relative or
%   absolute path may also be included in filename.
%
%   === Outputs ===
%   image_out is a 2D array of data from the specified file.
%
%   === Notes ===
%   This function is very low-level and may need to be adjusted for
%   different file formats.  As is, it is designed to read in a file of
%   tab-delimited data with no headers.  If formated this way, it can
%   handle files containing a rectangular array of arbitrary size.
%
%   === Example Usage ===
%   This is a helper function that cannot be called directly.

delimiter = '\t';

% Open the file.
filename=strtrim(filename); %Vendeiro trim whitespace
fileID = fopen(filename,'r');

%Read one line so we can figure out number of columns (knowing this
%increases the speed of textscan)
first_line=fgetl(fileID);

%Number of elements in a row is one more than number of delimiters
n_cols= sum(first_line(:)==char(sprintf(delimiter)))+1;

%Create formatSpec for interpretting rows
formatSpec = repmat('%f',[1,n_cols]); %makes '%f%f%f...'

%Rewind back to start of file before we read it all in
frewind(fileID);

%Read in the data with textscan.  It will return it as a cell array of cell
%arrays, with each sub-array corresponding to one column of data
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);

%Close the file
fclose(fileID);

%Convert data into a 2D array
image_out = [dataArray{1:end}];
end