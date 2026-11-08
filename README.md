# Absolutely continuous laws on a bounded interval are completely mixable for all large n

Yitzchak Shmalo

Paper: [paper/paper.pdf](paper/paper.pdf) (source [paper/paper.tex](paper/paper.tex) with
[paper/references.bib](paper/references.bib); amsart, pdflatex and bibtex).

A law $F$ is $n$-completely mixable if there are $n$ random variables with law $F$ whose sum is almost surely
constant. Puccetti, Wang and Wang proved this for laws on $[a,b]$ whose density is at least $3/(n(b-a))$ on the whole
interval, and Wang (Current open questions in complete mixability, Probab. Surveys 12 (2015), Open Problem 7) asked
whether every absolutely continuous law on a bounded interval is $n$-completely mixable for all large $n$.

The paper answers this in the affirmative, and more generally for every law on a bounded interval with a nonzero
absolutely continuous part. The main step is quantitative: a law on an interval of length $R$ that dominates $q$
times the uniform law on an interval of length $2H$ is $n$-completely mixable for every $n\ge(6R/H+1)/q$, whatever
the rest of the law is. The copies drawn from the rest of the law are coupled on a quantile grid, so that their sum
is within $R$ of its mean, and the error is absorbed by linear tilts of the uniform copies, which average back to the
uniform law. The general case groups the copies in pairs, since the law of the sum of two independent copies
dominates a measure with a continuous nonzero density; one extra copy handles odd $n$.

## Code

Plain Python 3 (standard library only).

- `code/check_general.py`: checks in exact rational arithmetic the identities used in the proof of Proposition 3.1
  (the quantile-grid bound and the mean of the grid sum, $|C|\le2R+L$, $|\beta|\le1$, $N-K\ge3$, $\mathbb E C=0$)
  for 1000 random laws with a discrete remainder, four widths of the extra law and three values of $N$. Output of
  the run: `code/check_general_run.txt`.
