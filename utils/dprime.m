function d = dprime(hitRate, falseAlarmRate)

d = norminv(hitRate)-norminv(falseAlarmRate);

return