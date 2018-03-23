%% This function generates a profiles with a random selection of attribute
%% levels varied from the BYO levels 

function [viewedProfiles,surveyquestions]=profile(viewedProfiles,min,max,numatt,resp,surveyAttrib,activeConstraints)

member=1;
  while member==1  %%Loop continues until a new profile is generated
  for i=1:numatt  %%
    surveyquestions(i,1)=resp.BYO(i,1);  %%profile contains all BYO
    vary(i,1)=round(min+rand*(max-min)); %%number of attributes to vary from BYO
    y=randsample(numatt,vary(i,1));  %%identify which attributes are varied
  for j=1:length(y)
    a=1:surveyAttrib(y(j,1),1).numberLevels;
    b=resp.BYO(y(j,1),1);
    if activeConstraints==1
    c=resp.revealedUnaccept(y(j,1),1);  %% revealed unacceptable are enforced
    else
    c=0;
    end
   surveyquestions(y(j,1),1)=randsample(setdiff(a,union(b,c)),1);  %%level of attribute
  end
      end
  if activeConstraints==1  %% must have are enforced
    for i=1:numatt  %% add back must have
      if resp.revealedMusthave(i,1)>0
      surveyquestions(i,1)=resp.musthave(i,1);
      end     
    end
  end
  member = ismember(transpose(surveyquestions(:,1)),transpose(viewedProfiles),'rows');
  end
 viewedProfiles=[viewedProfiles surveyquestions(:,1)];
end
  
  


