function [channel2, result]  =  sshfrommatlabissue(channel2,command)
%SSHFROMMATLAB issues commands to a remote computer from within Matlab
%
% [CONN, RESULT]  =  SSHFROMMATLABISSUE(CONN,COMMAND)
%
% Inputs:
%   CHANNEL is a Java ChannelShell object
%   COMMAND
% 
% Outputs:
%   CHANNEL is the returned Java ChannelShell object
%   RESULT
%
% See also SSHFROMMATLABCLOSE, SSHFROMMATLABINSTALL, SSHFROMMATLABISSUE
%
% (c) 2008 British Oceanographic Data Centre
%    Adam Leadbetter (alead@bodc.ac.uk)
%     2010 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 1.1
%
  import java.io.BufferedReader;
  import java.io.IOException;
  import java.io.InputStream;
  import java.io.InputStreamReader;
  import ch.ethz.ssh2.Connection;
  import ch.ethz.ssh2.Session;
  import ch.ethz.ssh2.StreamGobbler;
%
% Invocation checking
%
  if(nargin  ~=  2)
    error('Error: SSHFROMMATLABISSUE requires two input arguments...');
  end
  if(~isa(channel2,'ch.ethz.ssh2.Connection'))
    error(['Error: SSSHFROMMATLABISSUE input argument CHANNEL '...
      'is not a Java Connection object...']);
  end
  if(~ischar(command))
    error(['Error: SSSHFROMMATLABISSUE input argument COMMAND '...
      'is not a string...']);
  end
% 
% Send the commands
%
  result  =  {''};
  channel  =  channel2.openSession();
  channel.execCommand(command);
%
% Report the result to screen and to the string result...
%  
  stdout = StreamGobbler(channel.getStdout());
  br = BufferedReader(InputStreamReader(stdout));
  while(true)
    line = br.readLine();
	if(isempty(line))
	  break
    else
      if(isempty(result{1}))
        result{1}  =  char(line);
      else
        result{end+1}  =  char(line);
      end
	  %fprintf(1,'\n%s',char(line));
    end
  end
  channel.close();
  result  =  result';