let filename_showplace = document.getElementById("showplace");
filename_showplace.classList.add("max-h-60", "overflow-y-auto");
let stepped = false;
// check if the step is proccessed
document.addEventListener("DOMContentLoaded", async function () {
  const grabstep_rslt = await grabsteps(`/getSteps/create_adata/file_upload/`);
  console.log(grabstep_rslt);
  if (grabstep_rslt.message != "notfound") {
    stepped = true;
    document.getElementById("watchonly").classList.remove("hidden");
    filename_showplace.innerHTML = "";
    show_result_frontend(grabstep_rslt.output_values);
    for (let i = 0; i < grabstep_rslt.input_values.file_names.length; i++) {
      const file = {
        name: grabstep_rslt.input_values.file_names[i],
        size: grabstep_rslt.input_values.file_sizes[i],
      };
      show_file_frontend(file);
    }
    banned_operations();
  }
});
// can't operate
function banned_operations() {
  const processBtn = document.getElementById("processbutton");
  const fileUploadLabel = document.getElementById("file_upload_l");
  const delBtns = document.querySelectorAll(".delbtn");

  processBtn.disabled = true;
  fileUploadLabel.disabled = true;

  fileUploadLabel.classList.remove("bg-sky-700", "hover:bg-sky-600");
  fileUploadLabel.classList.add("bg-sky-700/50", "cursor-not-allowed");

  delBtns.forEach((btn) => (btn.disabled = true));
}
// proccess button click
async function processbtn(event) {
  event.preventDefault();
  const file_upload = document.getElementById("file_upload");
  const files = file_upload.files;

  toggleLoading(true, "processbutton");
  document.getElementById("adata_results").innerHTML = "";
  const errorElement = document.getElementById("error");
  const errorMessageElement = document.getElementById("errormessage");

  const formData = new FormData();

  for (let i = 0; i < files.length; i++) {
    formData.append("files", files[i]);
  }
  const csrftoken = getCookie("csrftoken");
  try {
    const result = await fetchAPI("/api/upload", formData, csrftoken);
    if (result.error) {
      throw new Error(result.error.message);
    }
    errorElement.classList.add("hidden");
    toggleLoading(false, "processbutton");
    show_result_frontend(result.data);
  } catch (error) {
    toggleLoading(false, "processbutton");
    nextbtn.classList.add("hidden");
    errorMessageElement.textContent = error.message || "上傳過程中發生錯誤";
    errorElement.classList.remove("hidden");
  }
}
// show result frontend
function show_result_frontend(result) {
  const ul = document.createElement("ul");
  ul.classList.add(
    "max-w-md",
    "space-y-1",
    "text-gray-500",
    "list-disc",
    "list-inside"
  );

  result.adata_results.forEach((item) => {
    const li = document.createElement("li");
    li.textContent = item;
    ul.appendChild(li);
  });

  adata_results.appendChild(ul);
  nextbtn.classList.remove("hidden");
}
// file upload
document.getElementById("file_upload").addEventListener("change", function (e) {
  filename_showplace.innerHTML = "";
  files = e.target.files;
  if (files.length > 0) {
    filename_showplace.classList.remove("justify-center");
    for (let i = 0; i < files.length; i++) {
      const file = files[i];
      show_file_frontend(file);
    }
  } else {
    filename_showplace.classList.add("justify-center");
    filename_showplace.textContent = "No Data Upload Yet.";
  }
});

// show file frontend
function show_file_frontend(file) {
  const filediv = document.createElement("div");
  filediv.classList.add(
    "flex",
    "w-full",
    "justify-between",
    "items-center",
    "mb-2",
    "p-2",
    "bg-gray-100",
    "rounded"
  );

  const filename = document.createElement("div");
  filename.textContent = file.name;
  filename.classList.add("truncate", "w-1/2");

  const filesize = document.createElement("div");
  filesize.textContent = formatFileSize(file.size);
  filesize.classList.add("text-sm", "text-gray-500");

  const delfile = document.createElement("button");
  delfile.textContent = "X";
  delfile.classList.add(
    "delbtn",
    "bg-red-500",
    "text-white",
    "px-2",
    "py-1",
    "rounded",
    "hover:bg-red-600",
    "disabled:bg-red-600/50",
    "disabled:cursor-not-allowed"
  );
  delfile.onclick = function () {
    filediv.remove();
    const file_upload = document.getElementById("file_upload");
    let files = Array.from(file_upload.files);

    let filename = file.name;

    files = files.filter((file) => file.name !== filename);
    file_upload.value = "";
    const dt = new DataTransfer();
    files.forEach((file) => dt.items.add(file));
    file_upload.files = dt.files;

    if (filename_showplace.children.length === 0) {
      filename_showplace.textContent = "No Data Upload Yet.";
    }
  };

  filediv.appendChild(filename);
  filediv.appendChild(filesize);
  filediv.appendChild(delfile);
  filename_showplace.appendChild(filediv);
}

function formatFileSize(bytes) {
  if (bytes < 1024) return bytes + " bytes";
  else if (bytes < 1048576) return (bytes / 1024).toFixed(2) + " KB";
  else if (bytes < 1073741824) return (bytes / 1048576).toFixed(2) + " MB";
  else return (bytes / 1073741824).toFixed(2) + " GB";
}
