function status = odeProgressReporter(t,~,~,reportTs)
if ~isempty(t)
    t=t(1);% for the intial function call
    if any(t==reportTs)
        runTimeSoFar=toc;
        runTimeSoFar=runTimeSoFar/60;
        disp([num2str(100*t/reportTs(end),'%.2f') '% of time points computed. runtime = ' num2str(runTimeSoFar,'%.2f') ' minutes'])
    end
end
status=0;
end