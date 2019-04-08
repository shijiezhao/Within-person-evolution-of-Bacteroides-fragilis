function send_jobs_to_cluster(cmds, jobsubmitoptions, Parallel, dirs)

%% Tami Lieberman, with contributions from Roy Kishony and Seungsoo Kim

% This  that runs commands on the unix shell (given in a cell array cmds)
% in parallel by generating temporary bash files
% that also produce output files and waitng for all output files to
% appear. Uses queue given in qname. Can run on mulitple dirs in a ingle
% command, Or multiple commands signle directory.

%Comes in 3 flavors, designated by Parallel
%Parallel = 1 -> run on cluster, wait for output
%Parallel = 2 -> run on cluster, do not wait for output
%Otherwise-- run sequentially, wait for output

%%

if isempty(cmds) | isempty(dirs)
    return
end



if Parallel==1
    dr = 'run_parallel_unix_commands_fast_tmp' ;
    outs={}; %list of orchestra output files to expect.
    if exist(dr,'dir')
        eval(['! rm ' dr '/*']) ;
    else
        mkdir(dr) ;
    end
    for i=1:max(length(cmds),length(dirs))
        fname = sprintf('%s/tmp%g.sh',dr,i) ;
        oname = sprintf('%s/out%g.txt',dr,i) ;
        outs{end+1}=oname;
        fid = fopen(fname,'w') ;
        fprintf(fid,'cd "%s"\n',dirs{min(i,end)}) ;
        fprintf(fid,'%s\n',cmds{min(i,end)}) ;
        fclose(fid) ;
        eval(sprintf('!chmod +x %s',fname))
       % eval(sprintf('!bsub -q %s ./%s',jobsubmitoptions,fname))
        eval(sprintf('!bsub -q %s -o %s ./%s',jobsubmitoptions,oname,fname))
        pause(1) % give the cluster a break
    end
    
    done = 0 ;
    while done < max(length(cmds),length(dirs))
        pause(10) ;
        for i=1:length(outs)
            if exist(outs{i},'file')
                fid=fopen(outs{i});
                l=fgetl(fid);
                while isempty(strfind(l,'Subject: Job '))
                    l=fgetl(fid);
                end
                if strfind(l,'Done')
                    done=done+1;
                    outs{i}=[];
                else
                    error(['A job failed. Check run_parallel_unix_commands_fast_tmp/out' num2str(i) '.txt for error message'])
                end
            end    
        end
        disp(done);
    end
    
elseif Parallel==2
    
    for i=1:max(length(cmds),length(dirs))
        fname = ['tmp' num2str(i) '.sh'] ;
        oname = ['out' num2str(i) '.txt'] ;
        fid = fopen(fname,'w') ;
        fprintf(fid,'cd "%s"\n',dirs{min(i,end)}) ;
        fprintf(fid,'%s\n',cmds{min(i,end)}) ;
        fclose(fid) ;
        eval(sprintf('!chmod +x %s',fname))
        eval(sprintf('!bsub -q %s -o %s ./%s',jobsubmitoptions,oname,fname))
    end
    fprintf('Continuing without waiting for last batch of jobs to finish...\n')
    
    
else
    
    cdr = pwd ;
    for i=1:max(length(cmds),length(dirs))
        fprintf(1, ['On sample ' num2str(i) ' ...\n']);
        cd(dirs{min(i,end)})
        eval(['!' cmds{min(i,end)}]) ;
        cd(cdr)
    end
end


end

