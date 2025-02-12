
//creates an HTTP request
function getHTTPObject()
{
 var xmlhttp;
 /*@cc_on
 @if (@_jscript_version >= 5)
  try {
   xmlhttp = new ActiveXObject("Msxml2.XMLHTTP");
    } catch (e) {
    try {
     xmlhttp = new ActiveXObject("Microsoft.XMLHTTP");
     } catch (E) {
       xmlhttp = false;
      }
    }
 @else
 xmlhttp = false;
 @end @*/
 if (!xmlhttp && typeof XMLHttpRequest != 'undefined') {
  try {
   xmlhttp = new XMLHttpRequest();
  } catch (e) {
    xmlhttp = false;
    }
 }
return xmlhttp;
}
var name_pw_http = "";
var end_sess = "";
var linka = new Array(10);
linka[0]=0;
linka[1]=0;
linka[2]=0;
linka[3]=0;
linka[4]=0;
linka[5]=0;
linka[6]=0;
linka[7]=0;
linka[8]=0;
linka[9]=0;
linka[10]=0;
linka[11]=0;
linka[12]=0;
var page_array = new Array("home", "calc", "db", "bb", "ab", "solar", "unit", "links", "", "faq",
"contact","ap","spectra");
//mover and mout function control the class names of the main menu bar.
function mover(mitem) {
 var menuitem;
 var menu_str = "link" + mitem;
 menuitem = document.getElementById(menu_str);
 var menu_class = menuitem.className;
 if(menu_class != "l_over") {
  menuitem.className="l_over";
 }
 document.body.style.cursor = 'pointer';
} //end function mover

function mout(mitem) {
 var menuitem;
 var mitem_str = "link" + mitem;
 document.body.style.cursor = 'default';
 menuitem = document.getElementById(mitem_str);
 var pg = document.getElementById("mypage").value;
 if(page_array[mitem] != pg) {
 menuitem.className="linactive";
 }
} //end function mout
//makes all of the main menu bars 'inactive'
function clear_all(item_num) {
 var menuitem;
 if(item_num != 100) {
 menuitem = document.getElementById("link" + item_num);
 menuitem.className="lactive";
 }
 var menuitems;
 for(var i=0; i<12; i++) {
  if(item_num != i)  {
   linka[i]= 0;
   menuitems = document.getElementById("link" + i);
   menuitems.className="linactive";
  }
  else  {
   if(item_num != 100)  {
   linka[i]=1;
   }
  }
 }
} //end function clearall


//this function retrieves the cookie from the users computer and sets the login name & password
//on the login page if the user has a valid cookie.
function load_cookie(){
 var ck_name = "";
 var ck_pw = "";
 var cke = Get_Cookie("spec_calc");
 if(cke.length > 0) {
  var ca = cke.split('-');
  ck_name = ca[0];
  ck_pw = ca[1];
  if((ck_name != "none") && (ck_name!="lim_user")) {
   var log_name = document.getElementById("spec_login");
   log_name.value = ck_name;
   var log_pw = document.getElementById("spec_pw");
   log_pw.value = ck_pw;
  }
 }
} //end load_cookie function

//this function gets the cookie from the users computer.
function Get_Cookie( name ) { 
 var arg= name + "=";
 var alen = arg.length;
 var clen = document.cookie.length;
 var i = 0;
 while(i < clen) {
  var j = i + alen; 
  if (document.cookie.substring(i, j) == arg)  {
   return getCookieVal(j);
  }
  i = document.cookie.indexOf(" " , i ) + 1;
  if (i == 0)  {
    break;
  }
 }
 return "none- ";
}
//gets the value from the stored cookie.
function getCookieVal (offset) {
 var endstr = document.cookie.indexOf (";", offset);
 if (endstr == -1) {
  endstr = document.cookie.length;
 }
 return unescape(document.cookie.substring(offset, endstr));
} //end getCookieVal
 
 //this function takes a name, value and # of days until expiration and sets a cookie on the users computer
function SetCookie(name,value,days){
 var date = new Date();
 date.setTime(date.getTime()+(days*24*60*60*1000));
 var expires = "; expires="+date.toGMTString();
 document.cookie = name+"="+value+expires+"; path=/";
}

//this function calls a php page to check the login name & pw.
function check_name_pw(){
 var spec_login_name = document.getElementById("spec_login").value;
 var spec_login_pw = document.getElementById("spec_pw").value;
 var login_spec_url = "../spec_login.php";
 var login_params = "action=checkpw&spec_login_name=" + spec_login_name + "&spec_login_pw=" + spec_login_pw;
 var remote_addr = document.getElementById("calc_remote_addr").value;
 var remote_host = document.getElementById("calc_remote_host").value;
 var user_agent = document.getElementById("calc_user_agent").value;
 login_params+= "&my_remote_address=" + remote_addr + "&my_remote_host=" + remote_host + "&my_http_agent=" + user_agent;
 name_pw_http = getHTTPObject();
 name_pw_http.open("POST", login_spec_url, true);
 name_pw_http.setRequestHeader("Content-Type", "application/x-www-form-urlencoded");
 name_pw_http.setRequestHeader("Content-length", login_params.length);
 name_pw_http.setRequestHeader("Connection", "close");
 name_pw_http.onreadystatechange = login_result;
 name_pw_http.send(login_params);
}

//gets the result from a users login and if a valid user & pw then enables all functionality and changes users status.
function login_result(){
 if(name_pw_http.readyState == 4) {
  var spec_login_name = document.getElementById("spec_login").value;
  var spec_login_pw = document.getElementById("spec_pw").value;
  var cookie_value = "";
  var log_res = name_pw_http.responseText;
  var res_array = log_res.split(":");
  var valid_login = res_array[0];
  var user_name = "";
  var user_id = "";
  if(valid_login == "true") {
   var pln_days = res_array[1];
   var uid = res_array[2];
   var ustat = res_array[4];
   var uplan = res_array[3];
   document.getElementById("user_name").value = spec_login_name;
   document.getElementById("user_id").value = uid;
   document.getElementById("login_err_div").innerHTML = "";
   document.getElementById("user_plan").value = uplan;
   user_name = spec_login_name;
   user_id = uid;
   var wl = window.location.href.toString();
   var wel_str = "";
   document.getElementById("user_status").value = ustat;
   if(ustat == "limited_reg") {
    wel_str = "<center><br><strong>Welcome " +  user_name + "!</strong><br><br>";
    wel_str+="  <a href='../info/account.php' class='logout_out2' onmouseover=className='logout_over2'; onmouseout=className='logout_out2';>Subscription Inactive</a> ";
    wel_str+="<br><a href='../info/account.php' class='logout_out2' onmouseover=className='logout_over2'; onmouseout=className='logout_out2'; target='_self'>My account</a><br><a onClick='logout();' return false;'  class='logout_out2' onmouseover=className='logout_over2'; onmouseout=className='logout_out2';><br>Log Out&nbsp;</a></center>";
    wel_str+="<span id='login_err_div'></span><br>";
    document.getElementById("namediv").innerHTML = wel_str;
   }
   else {
    wel_str = "<center><br><br><strong>Welcome " + user_name + "!</strong><br><br>";
    if((user_name != "aim_stm") && (user_name != "sensor") && (user_name != "limb") && (user_name != "limbscatter") && (user_name.toLowerCase() != "agu") && (user_name.toLowerCase() != "hitran") && (user_name.toLowerCase() != "spie")) {
    wel_str+="<a href='../info/account.php' class='logout_out2' onmouseover=className='logout_over2'; onmouseout=className='logout_out2'; target='_self'>My account</a>";
    }
    wel_str+="<br><a onClick='logout();' return false;' class='logout_out2' onmouseover=className='logout_over2'; onmouseout=className='logout_out2';><br>Log Out</a>&nbsp;</center>";
    wel_str+="<span id='login_err_div'></span><br>";
    document.getElementById("namediv").innerHTML = wel_str;
    cookie_value = spec_login_name + "-" + spec_login_pw;
    SetCookie('spec_calc',cookie_value,365);
    //document.getElementById("sub_pw_div").innerHTML = "";
    if(wl.indexOf('spectralcalc.php') != -1) {
     document.getElementById("ticks_box").checked = true;
     document.getElementById("calc_annotate").checked = true;
     document.getElementById("calc_timestamp").checked = true;
     document.getElementById("error_div").innerHTML = "";
     document.getElementById("instr_width").disabled = false;
     get_custom_linelist("calc_gas_db");
    } //end if on calculator tool
    else if(wl.indexOf('db_intensity.php') != -1) {
     get_custom_linelist("db_use_table");
    }
    else if(wl.indexOf('db_position.php') != -1) {
     get_custom_linelist("use_table");
    }
    else if(wl.indexOf('db_data.php') != -1) {
     get_custom_linelist("use_table");
    }
    else if(wl.indexOf('plans.php') != -1) {
     window.location = "../info/plans.php";
    }
    else if(wl.indexOf('spectra.php') != -1) {
     window.location = "../spectra/spectra.php";
    } 
    else if(wl.indexOf('account.php') != -1) {
     window.location = "../info/account.php";
    } 
    else if(wl.indexOf('modify_atmosphere.php') != -1) {
     //do nothing
    } 
    else if(wl.indexOf('atmosphere.php') != -1) {
     window.location = "../atmosphere_browser/atmosphere.php";
    } 
    else if(wl.indexOf('paths.php') != -1) {
    document.getElementById("path_submit").disabled = false;
    document.getElementById("reset_all").disabled = false;
    document.getElementById("clear_graph").disabled = false;
    document.getElementById("error_div").innerHTML = "";
    }
   } //end else
  } //end if valid_login = true;
  else
  {
   document.getElementById("login_err_div").innerHTML = "";
//   var err_div = document.getElementById("log_res_div");
  var err_div = document.getElementById("login_err_div");
  err_div.innerHTML = "<FONT COLOR='yellow' size='-1'>Invalid name/password.</font><br>";
  } //end else
 } //end if ready state = 4
} //end function login_result


function get_user_status() {
 var us = document.getElementById("user_status").value;
 return us;
}
function get_user_name() {
 var un = document.getElementById("user_name").value;
 return un;
}

function get_user_id() {
 var uid = document.getElementById("user_id").value;
 return uid;
}

//logout of the session.
function logout() {
 var remote_addr = document.getElementById("calc_remote_addr").value;
 var remote_host = document.getElementById("calc_remote_host").value;
 var user_agent = document.getElementById("calc_user_agent").value;
 var out="r_addr=" + remote_addr + "&r_host=" + remote_host + "&uag=" + user_agent;
 end_sess = getHTTPObject();
 end_sess.open("POST", "../logout_sess.php", true);
 end_sess.setRequestHeader("Content-Type", "application/x-www-form-urlencoded");
 end_sess.setRequestHeader("Content-length", out.length);
 end_sess.setRequestHeader("Connection", "close");
 end_sess.onreadystatechange = logout_result;
 end_sess.send(out); 
}

function logout_result(){
 if(end_sess.readyState == 4) {
  window.location = "../info/about.php";
 }
}

function chg_menu(this_tool) {
document.getElementById(this_tool).className = "mm_act";
// var md_str = "link" + this_tool;
// var hr_str = md_str + "_href";
// document.getElementById(md_str).className="l_over";
// document.getElementById(hr_str).className="menu_clr";
}

//this function gets a list of custom linelist that the user has associated with
//their user name or any group they belong to.
function get_custom_linelist(which_menu) {
 var params = "sel_menu=" + which_menu;
 ll_http = getHTTPObject();
 ll_http.open("POST", "../get_ll.php", true);
 ll_http.setRequestHeader("Content-Type", "application/x-www-form-urlencoded");
 ll_http.setRequestHeader("Connection", "close");
 ll_http.setRequestHeader("Content-length", params.length);
 ll_http.onreadystatechange = ll_res;
 ll_http.send(params);
} //end function get_custom_linelist

//this function updates the drop down menus to include any custom linelists.
function ll_res() {
 if(ll_http.readyState == 4) {
  var ll_xml = ll_http.responseXML;  
  var ll_tag = ll_xml.getElementsByTagName("ll").item(0);
  var sel_menu = ll_xml.getElementsByTagName("mu").item(0).firstChild.data;
  var ll_tbl = "";
  var ll_nme = "";
  var db_menu = document.getElementById(sel_menu);
  var numLL = ll_tag.childNodes.length;
  var db_len = db_menu.length; 
  with (db_menu){
   for(var q=0; q < numLL; q++) {
    ll_tbl =  ll_tag.getElementsByTagName("linelist").item(q).getAttribute("table_name");
    ll_nme = ll_tag.getElementsByTagName("linelist").item(q).firstChild.data;
    options[db_len] = new Option(ll_nme, ll_tbl);
    db_len++; 
   }
  }  //end with
 }
} //end function ll_res

function open_notice() {
window.open("../notice.html", "noticewWin", "statusbar=0, location=0, menubar=0, scrollbars=1, resizable=1, width=400, height=200");
}
