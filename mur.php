#!/usr/bin/php
<?php
# split -l 4974107 -a3 -d cactus-FCA-AJU_filtered.maf mafx/cactus.maf
# ./mur.php cactus-FCA-AJU.maf aju-ANGSD-all_updated.vcf

$maf = @$_SERVER['argv'][1];
$vcf = @$_SERVER['argv'][2];
$out = $maf . "_NEW_.res.vcf";

function REV($str)
{
    $str = strtoupper($str);
    $str = str_replace(["A","T","G","C"], [1,2,3,4], $str);
    $str = str_replace([2,1,4,3], ["A","T","G","C"], $str);
    return strrev($str);
}

function SAVE($vcf_data, $tofile = 'out.vcf')
{
    $H = fopen($tofile, 'w+');
    foreach ($vcf_data as $e => $list)
    {
        foreach ($list as $l)
        {
            if (count($l) == 16) continue ;
            if (count($l) == 16) { $l[] = "*"; $l[] = "*"; }
            fputs($H, implode(" ", $l) . "\n");
        }
    }
    #echo "Сохранено!\n";
    fclose($H);
}

# ---------------------------------------------------------------------------- #
echo "#1\n";

$vcf_data = [];
$H = fopen($vcf, "r");
if (!$H) die("\n.vcf not found\n");
while (($e = fgets($H)) !== false) 
{
    if ($e[0] == "#") continue;
    $e = str_replace("\n", "", $e);
    $x = explode(' ', $e);
    if (@!$vcf_data[$x[0]]) $vcf_data[$x[0]] = [];
    $vcf_data[$x[0]][] = $x;
}
fclose($H);

# ---------------------------------------------------------------------------- #
echo "#2\n";

$maf_data = [];
$H = fopen($maf, "r");
if (!$H) die("\n.maf not found\n");

$l = 0;
$first = false;
while (($e = fgets($H)) !== false) 
{
    $e = str_replace("\n", "", $e);
    #if ($l % 10000 == 0) echo "\r" . (100 * $l/4974107) . "%            ";
    if ($l % 10000 == 0) echo (100 * $l/4974107) . "%\n";
    if ($l % 10000 == 0) SAVE($vcf_data, $out);

    $l++;
    if (@$e[0] != "s") 
    {
        $first = false;
        continue;
    }
    if ($first) 
    {
        $box = [explode(' ', $first), explode(' ', $e)];
        sort($box);

        if ($box[0][4] == "-")
        {
            # self.global_start = self.all_length - self.start - self.length
            $box[0][2] = $box[0][5] - $box[0][2] - $box[0][3];
        }
        
        if ($box[0][1][0] == $box[1][1][0]) continue;
        $scaffold = str_replace("AcinonyxJubatus.", "", $box[0][1]);
        if (!@$vcf_data[$scaffold]) continue;

        foreach ($vcf_data[$scaffold] as $i => $m)
        {
            if ($m[1] <= $box[0][2]) continue ;
            if ($m[1] >= $box[0][2] + $box[0][3]) break ;
            
            $XX = strtoupper($box[0][6]);
            $YY = $box[1][6];

            # Перестраховка если гэпы
            $offset = strlen($XX) - strlen(str_replace("-", "", $XX));
            if ($m[1] < $box[0][2] + $box[0][3] - $offset) continue ; 


            if ($box[0][4] == "-")
            {
                $XX = REV($XX);
                $YY = REV($YY);
            }

            $m[3] = strtoupper($m[3]);
            $point = $m[1] - $box[0][2] - 1;
            
            for ($i = 0; $i < strlen($XX); $i++) if ($XX[$i] == "-") { $YY[$i] = "!"; $XX[$i] = "!"; }
            $XX = str_replace("!", "", $XX);
            $YY = str_replace("!", "", $YY);

            if ($XX[$point] == $m[3])
            {
                $m[] = substr($XX, $point - 1, 3);
                $m[] = substr($YY, $point - 1, 3);
                $vcf_data[$scaffold][$i] = $m;
            }
            else
            {
                # Этот сегмент можно игнорить
                # В идеальном мире мы тут не окажемся
                echo "ОТСТОЙ\n";
                echo substr($XX, $point-1, 3), " <- ", $m[3], "\n"; 
                print_r($box);
                print_r($m);
                echo "\n";
                if (rand(1,50) == 5) exit;
            }

        }
        $first = false;
    }
    else
    {
        $first = $e;
    }
}

SAVE($vcf_data, $out);
