
BEGIN {
    FS = "\t"
    OFS = "\t"
}

/^@HD/ {
    print
    next
}

/^@SQ/ {
    print
    next
}

/^@RG/ {
    print
    next
}

/^@PG/ {
    id = ""
    n = 0
    delete fields

    for (i = 1; i <= NF; i++) {
        if ($i ~ /^ID:/)
            id = substr($i, 4)

        if ($i !~ /^PP:/)
            fields[++n] = $i
    }

    if (id && seen[id]++)
        next

    for (i = 1; i <= n; i++)
        printf "%s%s", fields[i], (i == n ? "" : OFS)

    print ""
    next
}

{
    print
}
