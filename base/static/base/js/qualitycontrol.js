const f_sampleSelect = document.getElementById("f_sampleSelect");

// proccess button click
async function processbtn(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");

  const process_result = await fetchAPI("/api/qualitycontrol", 0, csrftoken);
  // init modal
  initmodal();
  //   show btn
  const btnbox = document.getElementById("btnbox");
  btnbox.classList.remove("hidden");
  //   show image
  showimage(process_result);
  // next page btn
  document.getElementById("nextbtn").classList.remove("hidden");
}

function initmodal() {
  document.getElementById("filterform").reset();
  document.getElementById("previewimagebox").innerHTML = "";
  document.getElementById("confirmbtn").classList.add("hidden");
}
// filter method change
function showimage(result) {
  f_sampleSelect.textContent = "";
  const imgbox = document.getElementById("imgbox");
  imgbox.textContent = "";
  const ul = document.createElement("ul");
  ul.classList.add(
    "space-y-1",
    "text-gray-500",
    "list-disc",
    "list-inside",
    "w-fit"
  );
  for (let i = 0; i < result.adata_results.length; i++) {
    const li = document.createElement("li");
    li.textContent = result["adata_results"][i];
    const img = document.createElement("img");
    img.src = `/media/tempimage/origin/${
      result["save_image_names"][i]
    }?t=${Date.now()}`;
    img.classList.add("m-5", "size-3/4");
    li.appendChild(img);
    ul.appendChild(li);

    // selection add data
    var sample = result["adata_results"][i].substring(
      0,
      result["adata_results"][i].indexOf(":")
    );
    var f_sampleoption = document.createElement("option");
    f_sampleoption.value = sample;
    f_sampleoption.textContent = sample;
    f_sampleSelect.appendChild(f_sampleoption);
  }
  imgbox.appendChild(ul);
}
document
  .getElementById("filter_method")
  .addEventListener("change", function () {
    console.log(this.value);
    let upandlow = document.getElementById("upandlow");
    let lo_limit_label = document.getElementById("lo_limit_label");
    let up_limit_label = document.getElementById("up_limit_label");
    if (this.value === "manual") {
      upandlow.classList.remove("hidden");
      lo_limit_label.textContent = "Please enter lower limit:";
      up_limit_label.textContent = "Please enter upper limit:";
    } else if (this.value === "quantile") {
      upandlow.classList.remove("hidden");
      lo_limit_label.textContent = "Please enter quantile of lower limit:";
      up_limit_label.textContent = "Please enter quantile of upper limit:";
    } else {
      upandlow.classList.add("hidden");
    }
  });

// preview
async function preview() {
  const previewimagebox = document.getElementById("previewimagebox");
  previewimagebox.textContent = "";
  const form = document.getElementById("filterform");
  const formData = new FormData(form);
  const csrftoken = getCookie("csrftoken");
  const preview_result = await fetchAPI("/api/preview", formData, csrftoken);
  // show preview image

  const ul = document.createElement("ul");
  ul.classList.add(
    "space-y-1",
    "text-gray-500",
    "list-disc",
    "list-inside",
    "w-fit"
  );
  const li = document.createElement("li");
  li.textContent = preview_result["adata_result"];
  const img = document.createElement("img");
  img.src = `/media/tempimage/preview/${
    preview_result["save_image_names"]
  }?t=${Date.now()}`;
  img.classList.add("m-5", "size-3/4");
  li.appendChild(img);
  ul.appendChild(li);

  previewimagebox.appendChild(ul);
  //   show confirm btn
  document.getElementById("confirmbtn").classList.remove("hidden");
}
async function confirm() {
  const csrftoken = getCookie("csrftoken");
  const replaceimg_result = await fetchAPI("/api/replace", 0, csrftoken);
  console.log(replaceimg_result);

  showimage(replaceimg_result);
}
// before next page
async function presubmit(event) {
  const csrftoken = getCookie("csrftoken");
  event.preventDefault();
  confirm_result = await fetchAPI("/api/confirm", 0, csrftoken);
  console.log(confirm_result);
  window.location.href = event.target.href;
}
// fetch api
async function fetchAPI(url, formData, csrftoken) {
  const response = await fetch(url, {
    method: "POST",
    body: formData,
    headers: {
      "X-CSRFToken": csrftoken,
    },
  });

  if (!response.ok) {
    throw new Error("Network response was not ok");
  }

  return response.json();
}
