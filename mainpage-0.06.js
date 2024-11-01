/*!
 * Code licensed under the Apache License v2.0.
 * For details, see http://www.apache.org/licenses/LICENSE-2.0.
 */

// jQuery for page scrolling feature - requires jQuery Easing plugin
$(function () {
    $(".page-scroll a").bind("click", function (event) {
        var $anchor = $(this);
        $("html, body")
            .stop()
            .animate(
                {
                    scrollTop: $($anchor.attr("href")).offset().top
                },
                500,
                "easeInOutExpo"
            );
        event.preventDefault();
    });
});

//resize canvas in jquery
(function ($) {
    $.fn.extend({
        //Let the user resize the canvas to the size he/she wants
        resizeCanvas: function (w, h) {
            var c = $(this)[0];
            c.width = w;
            c.height = h;
        }
    });
})(jQuery);

// Floating label headings for the contact form
$(function () {
    $("body")
        .on("input propertychange", ".floating-label-form-group", function (e) {
            $(this).toggleClass(
                "floating-label-form-group-with-value",
                !!$(e.target).val()
            );
        })
        .on("focus", ".floating-label-form-group", function () {
            $(this).addClass("floating-label-form-group-with-focus");
        })
        .on("blur", ".floating-label-form-group", function () {
            $(this).removeClass("floating-label-form-group-with-focus");
        });
});

// Highlight the top nav as scrolling occurs
$("body").scrollspy({
    target: ".navbar-fixed-top"
});

// Closes the Responsive Menu on Menu Item Click
$(".navbar-collapse ul li a").click(function () {
    $(".navbar-toggle:visible").click();
});

// show navbar once scrolled
var mainNavBarStatus = false;
var drawCreatureStatus = true;
var myp5;
$(document).scroll(function () {
    var height = $(document).scrollTop();
    var threshold = 20; // 20 pixel scroll, show/hide navbar.

    if (height <= threshold && mainNavBarStatus === true) {
        mainNavBarStatus = false;
        //$("#mainNavBar").fadeOut(500);
    } else if (height > threshold && mainNavBarStatus === false) {
        mainNavBarStatus = true;
        //$("#mainNavBar").fadeIn(500);
    }

    if (height > ($(window).height() * 3) / 4) {
        drawCreatureStatus = false;
    } else {
        drawCreatureStatus = true;
    }
});

function getContactMessage() {
    /*
  var currentdate = new Date();
  var datetime = "Sent on " + currentdate.getDay() + "/"+currentdate.getMonth() + "/" + currentdate.getFullYear() + " @ " + currentdate.getHours() + ":" + currentdate.getMinutes() + ":" + currentdate.getSeconds();
*/
    var datetime = new Date().toISOString().substr(0, 19);
    var result = {
        name:
            $("#contact_name").val() +
            "<br/><p></p><i>Sent on " +
            datetime +
            "</i>",
        email: $("#contact_email").val(),
        message: $("#contact_message").val(),
        read: false,
        archived: false
    };
    return result;
}

// this section deals with messages sent:
var SERVER_URL = "http://otoro.net:5984/otoro_contact/";
$("#contact_send").click(function () {
    $("#contact_success").hide();
    $("#contact_failure").hide();
    $("#contact_incomplete").hide();
    var theMessage = getContactMessage();
    // if incomplete:
    if (
        theMessage.name.length === 0 ||
        theMessage.email.length === 0 ||
        theMessage.message.length === 0
    ) {
        $("#contact_incomplete").fadeIn("slow");
        return;
    } else {
        // send the message out
        $.ajax({
            type: "POST",
            data: JSON.stringify(theMessage),
            async: true,
            url: SERVER_URL,
            contentType: "application/json",
            success: function (data) {
                var message = JSON.parse(data); // maybe useful later on in the future.
                $("#contact_success").fadeIn("slow");
            },
            error: function (e) {
                $("#contact_failure").fadeIn("slow");
            }
        });
    }
});

var mobileMode = false;
// main function
$(document).ready(function () {
    $("#mainNavBar").hide();
    $("#contact_success").hide();
    $("#contact_failure").hide();
    $("#contact_incomplete").hide();
    $("body").show();
    $("#creatureArea").show();

    mainNavBarStatus = false;

    console.log("page loaded.");

    var md = new MobileDetect(window.navigator.userAgent);
    if (md.mobile()) {
        mobileMode = true;
        console.log(md.mobile());
    } else {
        console.log("not mobile");
    }

    if (mobileMode) {
        var minDim = Math.round(
            1.0 * Math.min($(window).width(), $(window).height())
        );
        console.log("mobile dimention for creatures demo: " + minDim);
        $("#topCanvas").resizeCanvas(minDim, minDim);
        //        $.getScript('pendulum.mobile-0.06.js');
    } else {
        $("#topCanvas").resizeCanvas($(window).width(), $(window).height());
        //        $.getScript('pendulum-1.00.js');
    }
    //    $.getScript('js/p5.min.js');

    CREATURE.Main(mobileMode);
});
