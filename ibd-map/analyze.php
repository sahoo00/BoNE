<?php

header('Access-Control-Allow-Origin: *');

if (array_key_exists("go", $_GET)) {
  if (strcmp($_GET["go"], "listgenes") == 0) {
    listGenes($_GET);
  }
  if (strcmp($_GET["go"], "tindex") == 0) {
    findTindex($_GET);
  }
  if (strcmp($_GET["go"], "ibd01") == 0) {
    printIBD01($_GET);
  }
  if (strcmp($_GET["go"], "ibdOutcome") == 0) {
    printIBDOutcome($_GET);
  }
  if (strcmp($_GET["go"], "ibdGender") == 0) {
    printIBDGender($_GET);
  }
  if (strcmp($_GET["go"], "ibdMouse") == 0) {
    printIBDMouse($_GET);
  }
  if (strcmp($_GET["go"], "uc-map") == 0) {
    printUCMap($_GET);
  }
  if (strcmp($_GET["go"], "cd-map") == 0) {
    printCDMap($_GET);
  }
  if (strcmp($_GET["go"], "colon-tissue") == 0) {
    printColonTissue($_GET);
  }
  if (strcmp($_GET["go"], "ibdAll") == 0) {
    printAll($_GET);
  }
}

function getList($file, $index) {
    $glist = array();
    if (($fp = fopen($file, "r")) === FALSE) {
      echo "Can't open file $file <br>";
      exit;
    }
    $line = fgets($fp);
    while (!feof($fp))
    {
      $line = fgets($fp);
      $line = chop($line, "\r\n");
      $items = explode("\t", $line);
      if (count($items) > $index && $items[$index] != "") {
        array_push($glist, $items[$index]);
      }
    }
    fclose($fp);
    return $glist;
}

function listGenes($params) {
  $file = "genes.txt";
  $glist = getList($file, 0);
  sort($glist);
  echo json_encode($glist);
}

function processInput($str) {
  $list = preg_split("/\W/", $str);
  $res = implode(" ", $list);
  return $res;
}

function getOutput($cmd) {
    if ( ($fh = popen($cmd, 'r')) === false )
      die("Open failed: ${php_errormsg}\n");
    $data = 0;
    $imgurl = "";
    while (!feof($fh))
    {
      if ($data == 0) {
        $line = fgets($fh);
        $line = chop($line, "\r\n");
        if ($line == "--data--") {
          $data = 1;
        }
      }
      else {
        $line = fgets($fh);
        $imgurl .= $line;
      }
    }
    pclose($fh);
    return $imgurl;
}

function findTindex($list) {
  if (array_key_exists("id", $list)) {
    $id = processInput($list["id"]);
    $cmd = "perl analyze.pl tindex $id";
    $imgurl = getOutput($cmd);
    $type = "json";
    if (array_key_exists("type", $list)) {
      $type = $list["type"];
    }
    if ($imgurl != "") {
      if ($type == "html" || $type == "htmlpdf") {
        $json_a = json_decode($imgurl, true);
        echo $json_a['Minimum'] . '<br/>';
        echo $json_a['Maximum'] . '<br/>';
        echo $json_a['Total Num'] . '<br/>';
        echo '
<div id="regeStr"></div>
<script>
var jsonVar = ' . $imgurl . ';
var jsonStr = JSON.stringify(jsonVar);
var f = { brace: 0 }; 
var regeStr = "";
regeStr = jsonStr.replace(/({|}[,]*|[^{}:]+:[^{}:,]*[,{]*)/g, function (m, p1)
{
    var rtnFn = function() {
            return "<div style=\"text-indent: " + (f["brace"] * 20) + "px;\">" +
p1 + "</div>";
        },
        rtnStr = 0;
    if (p1.lastIndexOf("{") === (p1.length - 1)) {
        rtnStr = rtnFn();
        f["brace"] += 1;
    } else if (p1.indexOf("}") === 0) {
        f["brace"] -= 1;
        rtnStr = rtnFn();
    } else {
        rtnStr = rtnFn();
    }
    return rtnStr;
});

document.getElementById("regeStr").innerHTML += regeStr;
</script>
';
      }
      if ($type == "json") {
        echo $imgurl;
      }
    }
    else {
      if ($type == "html" || $type == "htmlpdf") {
        echo '<h2> Error </h2>';
      }
      if ($type == "json") {
        echo '{"error":1}';
      }
    }
  }
}

function printIBD01($list) {
  if (array_key_exists("id", $list)) {
    $id = processInput($list["id"]);
    $type = "json";
    if (array_key_exists("type", $list)) {
      $type = $list["type"];
    }
    $cmd = "HOME=tmpdir /usr/bin/python analyze.py --cmd ibd-1 --id \"$id\"";
    if ($type == "htmlpdf") {
      $cmd = "HOME=tmpdir /usr/bin/python analyze.py --cmd ibd-1 --pdf --id \"$id\"";
    }
    $imgurl = getOutput($cmd);
    if ($imgurl != "") {
      if ($type == "html") {
        $json_a = json_decode($imgurl, true);
        echo $json_a['roc'] . '<br/>';
        echo $json_a['roc-auc'] . '<br/>';
        echo $json_a['accuracy'] . '<br/>';
        echo $json_a['fisher'] . '<br/>';
        echo '<img src="data:image/jpeg;base64,' . $json_a['bar'] . '"/><br/>';
        echo '<img src="data:image/jpeg;base64,' . $json_a['density'] . '"/><br/>';
      }
      if ($type == "htmlpdf") {
        $json_a = json_decode($imgurl, true);
        echo $json_a['roc'] . '<br/>';
        echo $json_a['roc-auc'] . '<br/>';
        echo $json_a['accuracy'] . '<br/>';
        echo $json_a['fisher'] . '<br/>';
        echo '<iframe width="800px" height="800px" src="data:application/pdf;base64,' . $json_a['density'] . '"></iframe><br/>';
      }
      if ($type == "json") {
        echo $imgurl;
      }
    }
    else {
      if ($type == "html" || $type == "htmlpdf") {
        echo '<h2> Error </h2>';
      }
      if ($type == "json") {
        echo '{"error":1}';
      }
    }
  }
}

function printIBDOutcome($list) {
  if (array_key_exists("id", $list)) {
    $id = processInput($list["id"]);
    $type = "json";
    if (array_key_exists("type", $list)) {
      $type = $list["type"];
    }
    $cmd = "HOME=tmpdir /usr/bin/python analyze.py --cmd ibd-outcome --id \"$id\"";
    if ($type == "htmlpdf") {
      $cmd = "HOME=tmpdir /usr/bin/python analyze.py --cmd ibd-outcome --pdf --id \"$id\"";
    }
    $imgurl = getOutput($cmd);
    if ($imgurl != "") {
      if ($type == "html") {
        $json_a = json_decode($imgurl, true);
        for ($i = 0; $i < count($json_a['roc']); $i++) {
          echo $json_a['roc'][$i] . '<br/>';
          echo $json_a['roc-auc'][$i] . '<br/>';
          echo $json_a['accuracy'][$i] . '<br/>';
          echo $json_a['fisher'][$i] . '<br/>';
          echo $json_a['source'][$i] . '<br/>';
          echo $json_a['title'][$i] . '<br/>';
          echo $json_a['dbid'][$i] . '<br/>';
          echo $json_a['key'][$i] . '<br/>';
          echo '<img src="data:image/jpeg;base64,'.$json_a['bar'][$i].'"/><br/>';
        }
      }
      if ($type == "htmlpdf") {
        $json_a = json_decode($imgurl, true);
        for ($i = 0; $i < count($json_a['roc']); $i++) {
          echo $json_a['roc'][$i] . '<br/>';
          echo $json_a['roc-auc'][$i] . '<br/>';
          echo $json_a['accuracy'][$i] . '<br/>';
          echo $json_a['fisher'][$i] . '<br/>';
          echo $json_a['source'][$i] . '<br/>';
          echo $json_a['title'][$i] . '<br/>';
          echo $json_a['dbid'][$i] . '<br/>';
          echo $json_a['key'][$i] . '<br/>';
        }
        echo '<iframe width="800px" height="800px" src="data:application/pdf;base64,' . $json_a['outcome'] . '"></iframe><br/>';
      }
      if ($type == "json") {
        echo $imgurl;
      }
    }
    else {
      if ($type == "html" || $type == "htmlpdf") {
        echo '<h2> Error </h2>';
      }
      if ($type == "json") {
        echo '{"error":1}';
      }
    }
  }
}

function printIBDGender($list) {
  if (array_key_exists("id", $list)) {
    $id = processInput($list["id"]);
    $type = "json";
    if (array_key_exists("type", $list)) {
      $type = $list["type"];
    }
    $cmd = "HOME=tmpdir /usr/bin/python analyze.py --cmd ibd-gender --id \"$id\"";
    if ($type == "htmlpdf") {
      $cmd = "HOME=tmpdir /usr/bin/python analyze.py --cmd ibd-gender --pdf --id \"$id\"";
    }
    $imgurl = getOutput($cmd);
    if ($imgurl != "") {
      if ($type == "html") {
        $json_a = json_decode($imgurl, true);
        echo $json_a['source'] . '<br/>';
        echo $json_a['title'] . '<br/>';
        echo $json_a['dbid'] . '<br/>';
        echo $json_a['key'] . '<br/>';
        echo '<img src="data:image/jpeg;base64,' . $json_a['bar'] . '"/><br/>';
      }
      if ($type == "htmlpdf") {
        $json_a = json_decode($imgurl, true);
        echo $json_a['source'] . '<br/>';
        echo $json_a['title'] . '<br/>';
        echo $json_a['dbid'] . '<br/>';
        echo $json_a['key'] . '<br/>';
        echo '<iframe width="800px" height="800px" src="data:application/pdf;base64,' . $json_a['bar'] . '"></iframe><br/>';
      }
      if ($type == "json") {
        echo $imgurl;
      }
    }
    else {
      if ($type == "html" || $type == "htmlpdf") {
        echo '<h2> Error </h2>';
      }
      if ($type == "json") {
        echo '{"error":1}';
      }
    }
  }
}

function printIBDMouse($list) {
  if (array_key_exists("id", $list)) {
    $id = processInput($list["id"]);
    $type = "json";
    if (array_key_exists("type", $list)) {
      $type = $list["type"];
    }
    $cmd = "HOME=tmpdir /usr/bin/python analyze.py --cmd ibd-mouse --id \"$id\"";
    if ($type == "htmlpdf") {
      $cmd = "HOME=tmpdir /usr/bin/python analyze.py --cmd ibd-mouse --pdf --id \"$id\"";
    }
    $imgurl = getOutput($cmd);
    if ($imgurl != "") {
      if ($type == "html") {
        $json_a = json_decode($imgurl, true);
        for ($i = 0; $i < count($json_a['dbid']); $i++) {
          echo $json_a['roc-auc'][$i] . '<br/>';
          echo $json_a['accuracy'][$i] . '<br/>';
          echo $json_a['fisher'][$i] . '<br/>';
          echo $json_a['source'][$i] . '<br/>';
          echo $json_a['title'][$i] . '<br/>';
          echo $json_a['dbid'][$i] . '<br/>';
          echo $json_a['key'][$i] . '<br/>';
          echo $json_a['mtype'][$i] . '<br/>';
          echo $json_a['model'][$i] . '<br/>';
        }
        echo '<img src="data:image/jpeg;base64,'.$json_a['bar'].'"/><br/>';
      }
      if ($type == "htmlpdf") {
        $json_a = json_decode($imgurl, true);
        for ($i = 0; $i < count($json_a['dbid']); $i++) {
          echo $json_a['roc-auc'][$i] . '<br/>';
          echo $json_a['accuracy'][$i] . '<br/>';
          echo $json_a['fisher'][$i] . '<br/>';
          echo $json_a['source'][$i] . '<br/>';
          echo $json_a['title'][$i] . '<br/>';
          echo $json_a['dbid'][$i] . '<br/>';
          echo $json_a['key'][$i] . '<br/>';
          echo $json_a['mtype'][$i] . '<br/>';
          echo $json_a['model'][$i] . '<br/>';
        }
        echo '<iframe width="800px" height="800px" src="data:application/pdf;base64,' . $json_a['bar'] . '"></iframe><br/>';
      }
      if ($type == "json") {
        echo $imgurl;
      }
    }
    else {
      if ($type == "html" || $type == "htmlpdf") {
        echo '<h2> Error </h2>';
      }
      if ($type == "json") {
        echo '{"error":1}';
      }
    }
  }
}

function printUCMap($list) {
  if (array_key_exists("id", $list)) {
    $id = processInput($list["id"]);
    $cmd = "perl analyze.pl uc-map $id";
    $obj = getOutput($cmd);
    $type = "json";
    if (array_key_exists("type", $list)) {
      $type = $list["type"];
    }
    if ($obj != "") {
      if ($type == "html" || $type == "htmlpdf") {
        $json_a = json_decode($obj, true);
        $ll = explode(" ", $id);
        for ($i = 0; $i < count($ll); $i++) {
          if (array_key_exists($ll[$i], $json_a)) {
            echo $json_a[$ll[$i]] . '<br/>';
          }
        }
      }
      if ($type == "json") {
        echo $obj;
      }
    }
    else {
      if ($type == "html" || $type == "htmlpdf") {
        echo '<h2> Error </h2>';
      }
      if ($type == "json") {
        echo '{"error":1}';
      }
    }
  }
}

function printCDMap($list) {
  if (array_key_exists("id", $list)) {
    $id = processInput($list["id"]);
    $cmd = "perl analyze.pl cd-map $id";
    $obj = getOutput($cmd);
    $type = "json";
    if (array_key_exists("type", $list)) {
      $type = $list["type"];
    }
    if ($obj != "") {
      if ($type == "html" || $type == "htmlpdf") {
        $json_a = json_decode($obj, true);
        $ll = explode(" ", $id);
        for ($i = 0; $i < count($ll); $i++) {
          if (array_key_exists($ll[$i], $json_a)) {
            echo $json_a[$ll[$i]] . '<br/>';
          }
        }
      }
      if ($type == "json") {
        echo $obj;
      }
    }
    else {
      if ($type == "html" || $type == "htmlpdf") {
        echo '<h2> Error </h2>';
      }
      if ($type == "json") {
        echo '{"error":1}';
      }
    }
  }
}

function printColonTissue($list) {
  if (array_key_exists("id", $list)) {
    $id = processInput($list["id"]);
    $cmd = "HOME=tmpdir /usr/bin/python analyze.py --cmd colon-tissue --id \"$id\"";
    $obj = getOutput($cmd);
    $type = "json";
    if (array_key_exists("type", $list)) {
      $type = $list["type"];
    }
    if ($obj != "") {
      if ($type == "html" || $type == "htmlpdf") {
        $json_a = json_decode($obj, true);
        $ll = explode(" ", $id);
        for ($i = 0; $i < count($ll); $i++) {
          if (array_key_exists($ll[$i], $json_a)) {
            echo $ll[$i]." ";
            print_r($json_a[$ll[$i]]);
            echo '<br/>';
          }
        }
      }
      if ($type == "json") {
        echo $obj;
      }
    }
    else {
      if ($type == "html" || $type == "htmlpdf") {
        echo '<h2> Error </h2>';
      }
      if ($type == "json") {
        echo '{"error":1}';
      }
    }
  }
}

function printAll($list) {
  echo '<h2> Therapeutic Index </h2>';
  findTindex($list);
  echo '<h2> Analysis of GSE83687 Peters 2017 IBD dataset </h2>';
  printIBD01($list);
  echo '<h2> Analysis of IBD Outcome dataset </h2>';
  printIBDOutcome($list);
  echo '<h2> Analysis of IBD Gender Bias </h2>';
  printIBDGender($list);
  echo '<h2> Analysis of IBD Mouse Models</h2>';
  printIBDMouse($list);
  echo '<h2> Analysis of UC map</h2>';
  printUCMap($list);
  echo '<h2> Analysis of CD map</h2>';
  printCDMap($list);
  echo '<h2> Analysis of Colon Tissue expression profile</h2>';
  printColonTissue($list);
}


?>
