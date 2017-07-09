function csvwriteFast( filename, values, format )
% values: numeric matrix
% format: C format, separate for each column of values.
% example: 
% r=randn(1e5,3);
% csvwriteFast( 'r.csv', r, '%f,%f,%f' )

format = [ format, '\n' ];
fid = fopen( filename, 'w' );
fprintf( fid, format, values' );
fclose( fid );
return
