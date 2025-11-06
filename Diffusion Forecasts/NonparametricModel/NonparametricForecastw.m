function [ferr,Forecast,Truth,ForecastVar] = NonparametricForecastw(trainingData,obs,truth,R,forecastSteps,delays,nvars)

%%% Written by Tyrus Berry of George Mason University 
%%% --> augmented for weighted methods by Jason J. Bramburger

    N = size(trainingData,1);
    M = size(trainingData,2);
    
    verSteps = size(obs,1);
    initstep = 1;
    verLength = (verSteps-initstep-forecastSteps);

    Forecast = zeros(verLength,forecastSteps,M);
    Truth = zeros(verLength,forecastSteps,M);
    ForecastVar = zeros(verLength,forecastSteps,M);
    ferr = zeros(forecastSteps,1);

    mdelay=max(delays);
    T = N-mdelay+1;

    %%%%%%%%%%%%%%%%%%%%%% Build Model %%%%%%%%%%%%%%%%%%%%%%

    k=512;
    k2=22;
    shifts=1;
    [basis,ForwardOp,~,peq,qest,normalizer] = BuildModelVBDMw(trainingData,k,k2,nvars,delays,shifts);
    xPred = trainingData(mdelay:end,:);

    %%%%%%%%%%%%%%%%%%%%%% Filter and Forecast %%%%%%%%%%%%%%%%%%%%%%
    
    if (R==0)
        recon = basis*((basis)'*(xPred.*repmat(peq./qest,1,size(xPred,2))));
        rerrs = recon - xPred;
        R = rerrs'*rerrs/size(rerrs,1);
    end    
    [u,s] = svd(R);
    invrootR = u*diag(1./sqrt(diag(s)))*u';

    meanCoeff = (repmat(peq./qest,1,M).*xPred)'*basis/T;
    varCoeff = (repmat(peq./qest,1,M).*(xPred.^2))'*basis/T;

    %%% Initialize on the invariant measure
    p0 = peq;

    %%% Project onto the eigenfunction basis
    c = (basis)'*(p0./qest)/T;
    c = c/(normalizer*c);

    for i = 1:verSteps

        %%% Forecast Step of filter
        c = ForwardOp*c;
        c = c/(normalizer*c);
 
        %%% Bayesian filter
        llhood = exp(-sum((invrootR*(repmat(obs(i,:),T,1)-xPred)').^2,1))';
        prior = (repmat(peq,1,nvars).*basis)*c;
        post = max((prior).*llhood,0);
        %post = max((prior+peq*fparam).*llhood,0);
        post = post./mean(post./qest);
        c = (basis)'*(post./qest)/T;
        c = c/(normalizer*c);
        
        %%% Forecast from the posterior estimate
        ctempj = c;
        if ((i>=initstep)&&(i<verSteps-forecastSteps-1))     
            for j=1:forecastSteps

                Forecast(i,j,:) = meanCoeff*ctempj;
                Truth(i,j,:) = truth(i+j-1,:);
                
                ferr(j) = ferr(j) + mean((truth(i+j-1,:) - squeeze(Forecast(i,j,:))').^2,2)/verLength;
                
                ForecastVar(i,j,:) = sqrt(abs(varCoeff*ctempj - Forecast(i,j,:)'.^2));

                ctempj = ForwardOp*ctempj;
                ctempj = ctempj/(normalizer*ctempj); 

            end
        end
    end
    
    ferr = sqrt(ferr);

end

