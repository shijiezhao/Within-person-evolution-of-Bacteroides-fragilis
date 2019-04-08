function p=ChiSqPval(seq1,seq2,n);
expected = seq1/sum(seq1);
expected = expected*sum(seq2);
chi = sum((seq2-expected).^2./expected);
p= 1-chi2cdf(chi,n);
end