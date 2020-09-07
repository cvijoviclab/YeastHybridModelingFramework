function [path] = createResultsFolder(text, var, gene)
%written by: Linnea Österberg
%date: 2020-04-02
%description: With given input creates a directory for saving results in 
%the results folder.
%input: string with folder name (1) and variable defining this experiment
%(2)

foldername = ([text, (num2str(var)) '/' char(gene)]);
       
%% create folders if they don't already exist
if isfolder( 'results')== 0
    mkdir results;
end

i=2;
if isfolder(['results/',foldername]) == 0
    mkdir('results', foldername);
    path = ['results/', foldername,'/'];
else 
    newFoldername=foldername;
    while isfolder(['results/',newFoldername]) == 1
        newFoldername = [foldername,'nr',num2str(i)];
        i=i+1;
    end
    mkdir('results', newFoldername);
    path = ['results/', newFoldername,'/'];
end


end