classdef respondent
 properties
    BYO
    partworth
    musthave
    unaccept
    revealedMusthave
    revealedUnaccept
    surveyquestions
    surveyresponses
    surveydata
    tournament
   end

 methods
 function decision=evaluate(obj, qlevels,numAttributes,numProfiles)
  for i=1:numProfiles
   decision(1,i)=1;
          for j=1:numAttributes
           if obj.musthave(j,1)~=0   %%check that musthave are observed
           if qlevels(j,i)==obj.musthave(j,1) 
            decision(1,i)=1*decision(1,i);
          else
         decision(1,i)=0;
          end
        end
      end
          for j=1:numAttributes
           if obj.unaccept(j,1)~=0    %%check that unacceptable are observed
             if qlevels(j,i)~=obj.unaccept(j,1)
             decision(1,i)=1*decision(1,i);
             else
             decision(1,i)=0;
             end
           end
       end
  end
 end
    
            
function obj=processChoices(obj,numAttributes,numProfiles,firstquestion,lastquestion)  
        
  for currentquestion=firstquestion:lastquestion
    for i=1:numProfiles
      for j=1:numAttributes
           if obj.surveydata((currentquestion-1)*numAttributes+j,2+numProfiles+i)==1
              if obj.musthave(j,1)==obj.surveydata((currentquestion-1)*numAttributes+j,i)
                  obj.revealedMusthave(j,1)=obj.musthave(j,1);
              end
           end
           if obj.surveydata((currentquestion-1)*numAttributes+j,2+numProfiles+i)==0
              if obj.unaccept(j,1)==obj.surveydata((currentquestion-1)*numAttributes+j,i)
               obj.revealedUnaccept(j,1)=obj.unaccept(j,1);
              end
           end
      end    
    end
   end
end
end
end




