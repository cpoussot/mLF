function latexList = make_latex_examples()

%
COUNT_POLY  = 0;
COUNT_RAT   = 0;
COUNT_IRRAT = 0;
%
latexList   = [];
%latexList   = '\small ';
latexList   = [latexList '\begin{center} \begin{longtable}{|c|l|l|} \hline '];
latexList   = [latexList 'Case & Info. & Function \\ \hline\hline '];
for CAS = 1:48
    COLOR       = '';
    [H,info]    = mlf.examples(CAS);
    for i = 1:numel(info.tag)
        if strcmp(info.tag{i},'polynomial')
            COUNT_POLY  = COUNT_POLY+1;
            COLOR       = '\blue';
        end
        if strcmp(info.tag{i},'rational')
            COUNT_RAT = COUNT_RAT+1;
        end
        if strcmp(info.tag{i},'irrational')
            COUNT_IRRAT = COUNT_IRRAT+1;
            COLOR       = '\orange';
        end
    end
    latexList = [latexList ...
                 COLOR '{\#' num2str(CAS) '} & ' ... 
                 ' $\ord= ' num2str(info.n) '$, tensor size: ' num2str(info.tab_MB,3) ' MB & ' ... 
                 info.name ' \\' ];
end
latexList = [latexList '\hline \caption{List of examples. \blue{Polynomial: ' num2str(COUNT_POLY) '}, {rational: ' num2str(COUNT_RAT) '}, \orange{irrational: ' num2str(COUNT_IRRAT) '}.} \label{tab:examples} \end{longtable} \end{center}\normalsize' ];
