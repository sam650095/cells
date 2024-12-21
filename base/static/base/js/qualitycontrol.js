const f_sampleSelect = document.getElementById("f_sampleSelect");
let marker_list = [];
// check if the step is proccessed
document.addEventListener("DOMContentLoaded", async function () {
  const grabstep_rslt = await grabsteps(`/getSteps/qualitycontrol/process/`);
  const grabstep_filter_rslt = await grabsteps(
    `/getSteps/qualitycontrol/filter/`
  );
  console.log(grabstep_rslt);
  console.log(grabstep_filter_rslt);
  if (
    grabstep_rslt.message != "notfound" ||
    grabstep_filter_rslt.message != "notfound"
  ) {
    stepped = true;
    document.getElementById("watchonly").classList.remove("hidden");
    document.getElementById("btnbox").classList.remove("hidden");
    document.getElementById("nextbtn").classList.remove("hidden");
    showimage(grabstep_rslt.output_values);
    if (grabstep_filter_rslt.message != "notfound") {
      document.getElementById("sampleul_select").textContent =
        grabstep_filter_rslt.input_values["f_sampleSelect"];
      document.getElementById("f_sampleSelect").disabled = true;
      document.getElementById("sampleul_input").value =
        grabstep_filter_rslt.input_values["f_sampleSelect"];
      document.getElementById("sampleul_input").disabled = true;
      document.getElementById("minGenes").value =
        grabstep_filter_rslt.input_values["minGenes"];
      document.getElementById("minGenes").disabled = true;
      document.getElementById("filter_sampleul_select").textContent =
        grabstep_filter_rslt.input_values["filter_method"];
      document.getElementById("filter_method").disabled = true;
      document.getElementById("filter_sampleul_input").value =
        grabstep_filter_rslt.input_values["filter_method"];
      document.getElementById("filter_sampleul_input").disabled = true;
      select_method(grabstep_filter_rslt.input_values["filter_method"]);
      document.getElementById("lowerlimit").value =
        grabstep_filter_rslt.input_values["lowerlimit"];
      document.getElementById("lowerlimit").disabled = true;
      document.getElementById("upperlimit").value =
        grabstep_filter_rslt.input_values["upperlimit"];
      document.getElementById("upperlimit").disabled = true;
      document.getElementById("confirmbtn").classList.remove("hidden");
      if (grabstep_filter_rslt.output_values["adata_result"] != undefined) {
        show_preview_image(grabstep_filter_rslt.output_values);
      }
      banned_operations("filter");
    }
    banned_operations("process");
  }
});
// can't operate
function banned_operations(m) {
  const processBtn = document.getElementById("processbutton");
  processBtn.disabled = true;
  const previewBtn = document.getElementById("previewbutton");
  const confirmBtn = document.getElementById("confirmbtn");
  processBtn.disabled = true;
  previewBtn.disabled = true;
  confirmBtn.disabled = true;
}
// proccess button click
async function processbtn(event) {
  event.preventDefault();
  toggleLoading(true, "processbutton");
  const csrftoken = getCookie("csrftoken");
  try {
    const process_result = await fetchAPI("/api/qualitycontrol", 0, csrftoken);
    if (process_result.error) {
      throw new Error(process_result.error.message);
    }
    // init modal
    initmodal();
    //   show btn
    document.getElementById("btnbox").classList.remove("hidden");
    //   show image
    showimage(process_result.data);
    console.log(process_result);
    toggleLoading(false, "processbutton");
    // next page btn
    document.getElementById("nextbtn").classList.remove("hidden");
    marker_list = process_result.data.marker_list;
    console.log(marker_list);
  } catch (error) {
    toggleLoading(false, "processbutton");
  }
}

function initmodal() {
  document.getElementById("filterform").reset();
  document.getElementById("previewimagebox").innerHTML = "";
  document.getElementById("confirmbtn").classList.add("hidden");
}
const toggleDropdown = (e) => {
  e.preventDefault();
  setIsOpen(!isOpen);
};
// filter method change
function showimage(result) {
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

  const sample_ul = document.getElementById("sampleul");
  sample_ul.textContent = "";
  defaultsample = result["adata_results"][0].substring(
    0,
    result["adata_results"][0].indexOf(":")
  );
  document.getElementById("sampleul_select").textContent = defaultsample;
  document.getElementById("sampleul_input").value = defaultsample;
  for (let i = 0; i < result.adata_results.length; i++) {
    const li = document.createElement("li");
    li.textContent = result["adata_results"][i];
    const img = document.createElement("img");
    img.src = `/media/qualitycontrol/origin/${
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

    const sample_li = document.createElement("li");
    sample_li.innerHTML = `
              <a onclick="select('sampleul','${sample}')"
              class="block px-4 py-2 hover:bg-gray-100 dark:hover:bg-gray-600 dark:hover:text-white m-2 cursor-pointer"
          >${sample}</a>
      `;

    sample_ul.appendChild(sample_li);
  }
  imgbox.appendChild(ul);
}
// change dropdown for selecting filter method
function select_method(method) {
  select("filter_sampleul", method);
  let upandlow = document.getElementById("upandlow");
  let lo_limit_label = document.getElementById("lo_limit_label");
  let lowerlimit = document.getElementById("lowerlimit");
  let up_limit_label = document.getElementById("up_limit_label");
  let upperlimit = document.getElementById("upperlimit");

  if (method == "Manual") {
    upandlow.classList.remove("hidden");
    lo_limit_label.textContent = "Please enter lower limit:";
    up_limit_label.textContent = "Please enter upper limit:";
    lowerlimit.value = "";
    upperlimit.value = "";
  } else if (method == "Quantile") {
    upandlow.classList.remove("hidden");
    lo_limit_label.textContent = "Please enter quantile of lower limit:";
    up_limit_label.textContent = "Please enter quantile of upper limit:";
    lowerlimit.value = "";
    upperlimit.value = "";
  } else {
    upandlow.classList.add("hidden");
    lowerlimit.value = "";
    upperlimit.value = "";
  }
}
function select(s_ul, selected) {
  document.getElementById(s_ul + "_select").textContent = selected;
  document.getElementById(s_ul + "_input").value = selected;
}
// preview
async function preview() {
  const previewimagebox = document.getElementById("previewimagebox");
  previewimagebox.textContent = "";
  const form = document.getElementById("filterform");
  const formData = new FormData(form);
  const csrftoken = getCookie("csrftoken");
  const preview_result = await fetchAPI("/api/preview", formData, csrftoken);
  console.log(preview_result);
  show_preview_image(preview_result.data);
}
// show preview image
function show_preview_image(preview_result) {
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
  img.src = `/media/qualitycontrol/preview/${
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
  const form = document.getElementById("filterform");
  const formData = new FormData(form);
  const replaceimg_result = await fetchAPI("/api/replace", formData, csrftoken);

  showimage(replaceimg_result.data);
}
// before next page
async function presubmit(event) {
  const csrftoken = getCookie("csrftoken");
  event.preventDefault();
  confirm_result = await fetchAPI("/api/confirm", 0, csrftoken);
  window.location.href = event.target.href;
}
